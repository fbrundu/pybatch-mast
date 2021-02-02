"""Main module."""

from typing import Optional, Generator, Tuple, List, Sequence, Dict, Any, Union

from anndata import AnnData
from pandas import DataFrame

import boto3 as bt
from botocore.exceptions import ClientError
import os
import pandas as pd
import scanpy as sc
import tempfile
import time
import uuid

class MASTCollectionError(Exception):
    def __init__(
        self,
        message: str,
        job_collection,
    ):
        # Call the base class constructor with the parameters it needs
        super().__init__(message)
        self.jc = job_collection


class BatchMAST():
    def __init__(
        self,
        job_queue: str,
        job_def: str,
        bucket: str,
        layer: str = 'counts',
    ):
        self.job_queue = job_queue
        self.job_def = job_def
        self.bucket = bucket
        self.layer = layer

    def mast(
        self,
        adata: AnnData,
        keys: Sequence[str],
        group: str,
        fdr: float,
        lfc: float,
        covs: str = '',
        bys: Optional[Sequence[Tuple[str, Sequence[str]]]] = None,
        min_perc: Optional[Union[float, Dict[str, float]]] = None,
        on_total: Optional[bool] = False,
        min_cells_limit: Optional[int] = 3,
        jobs: int = 1,
    ) -> Generator[
        Tuple[
            Dict[str, DataFrame],
            Dict[str, Dict[str, List[str]]],
            Optional[str],
        ],
        None,
        None,
    ]:
        # NOTE n_genes is always assumed as covariate
        if bys is None:
            if min_perc is not None:
                adata = adata.copy()
                if on_total:
                    total_cells = adata.shape[0]
                else:
                    total_cells = adata.obs[group].value_counts().min()
                min_cells = max(total_cells * min_perc, min_cells_limit)
                print(
                    f'Filtering genes detected in fewer than {min_cells} cells'
                )
                sc.pp.filter_genes(adata, min_cells=min_cells)
            enough_genes = adata.shape[1] > 0
            job_collection = {}
            if enough_genes:
                job_collection = self._mast(
                    job_collection, adata, covs, group, keys, jobs=jobs,
                )
            else:
                print('Not enough genes, computation skipped')
            try:
                de, top = self.mast_prep_output(job_collection, lfc, fdr)
            except ClientError as e:
                raise MASTCollectionError(e, job_collection) from e
            except Exception as e:
                raise MASTCollectionError(e.message, job_collection) from e
            yield de, top, None
        else:
            for by, groups in bys:
                job_collection = {}
                for b in groups:
                    adata_b = adata[adata.obs[by] == b].copy()
                    if min_perc is not None:
                        if on_total:
                            total_cells = adata_b.shape[0]
                        else:
                            total_cells = adata_b.obs[group].value_counts(
                            ).min()
                        min_cells = max(
                            total_cells * min_perc[b], min_cells_limit
                        )
                        print(
                            'Filtering genes detected in fewer '
                            f'than {min_cells} cells'
                        )
                        sc.pp.filter_genes(
                            adata_b, min_cells=min_cells,
                        )
                    enough_groups = (
                        adata_b.obs[group].value_counts() >= 3
                    ).sum() > 1
                    enough_genes = adata_b.shape[1] > 0
                    if enough_groups and enough_genes:
                        job_collection = self._mast(
                            job_collection, adata_b, covs, group,
                            keys, by=by, b=b, jobs=jobs,
                        )
                    else:
                        print(f'Computation for {b} skipped')
                try:
                    de, top = self.mast_prep_output(job_collection, lfc, fdr)
                except ClientError as e:
                    raise MASTCollectionError(e, job_collection) from e
                except Exception as e:
                    raise MASTCollectionError(e.message, job_collection) from e
                yield de, top, by

    def _mast(
        self,
        job_collection: Dict[str, Dict[str, str]],
        adata: AnnData,
        covs: str,
        group: str,
        keys: Sequence[str],
        by: Optional[str] = None,
        b: Optional[str] = None,
        jobs: int = 1,
    ) -> Dict[str, Dict[str, str]]:
        if by is None:
            b = 'Sheet0'
        new_covs = BatchMAST._clean_covs(adata, covs, group, by=by)
        remote_dir, job_id, job_name, content = self.mast_compute(
            adata, keys, group=group, covs=new_covs, block=False, jobs=jobs,
        )
        job_collection[job_id] = {'group': b, 'remote_dir': remote_dir}
        return job_collection

    def mast_prep_output(
        self,
        job_collection: Dict[str, Dict[str, str]],
        lfc: float,
        fdr: float,
        wait: float = 30,
    ) -> Tuple[DataFrame, Dict[str, Dict[str, List[str]]]]:
        de = {}
        top = {}
        for job_id, status, metadata, content in self.mast_collect(
            job_collection, wait=wait,
        ):
            if status == 'SUCCEEDED':
                b = metadata['group']
                de[b] = content
            elif status == 'FAILED':
                print(f'Job Failed: group {metadata["group"]}')
            else:
                raise NotImplementedError(f'Status {status} not managed')
        top = BatchMAST.mast_filter(de, lfc, fdr)
        return de, top

    @staticmethod
    def mast_filter(
        de: Dict[str, DataFrame],
        lfc: float,
        fdr: float,
    ) -> Dict[str, Dict[str, List[str]]]:
        top = {}
        for b in de.keys():
            cols = [
                '_'.join(c.split('_')[:-1])
                for c in de[b].columns[de[b].columns.str.endswith('_coef')]
            ]
            top[b] = {}
            for c in cols:
                top[b][c] = de[b][
                    (de[b][f'{c}_fdr'] < fdr) & (de[b][f'{c}_coef'] > lfc)
                ].sort_values([
                    f'{c}_fdr', f'{c}_coef'
                ], ascending=[True, False]).index.tolist()
        return top

    @staticmethod
    def _clean_covs(
        adata: AnnData,
        covs: str,
        group: str,
        by: Optional[str] = None,
    ) -> str:
        covs_s = covs.split('+')[1:]
        new_covs = ''
        for c in covs_s:
            if c not in (group, by) and adata.obs[c].nunique() > 1:
                new_covs += f'+{c}'
        return new_covs

    def _mast_prep(
        self,
        adata: AnnData,
        remote_dir: str,
        keys: Sequence[str],
        group: str,
        covs: str = '',
        ready: Optional[Sequence[str]] = None,
        jobs: int = 1,
    ) -> str:
        s3 = bt.resource('s3')
        if ready is None:
            ready = []
        with tempfile.TemporaryDirectory() as td:
            if 'mat' not in ready:
                local_mat = os.path.join(td, 'mat.fth')
                adata = adata.copy()
                adata.X = adata.layers[self.layer]
                sc.pp.normalize_total(adata, target_sum=1e6)
                sc.pp.log1p(adata, base=2)
                adata.to_df().reset_index().to_feather(
                    local_mat, compression='uncompressed',
                )
                remote_mat = os.path.join(remote_dir, 'mat.fth')
                print(f'Uploading matrix ({adata.shape}) to s3...')
                s3.meta.client.upload_file(local_mat, self.bucket, remote_mat)

            if 'cdat' not in ready:
                local_cdat = os.path.join(td, 'cdat.csv')
                adata.obs[keys].to_csv(local_cdat)
                remote_cdat = os.path.join(remote_dir, 'cdat.csv')
                print('Uploading metadata to s3...')
                s3.meta.client.upload_file(
                    local_cdat, self.bucket, remote_cdat)

            remote = os.path.join(self.bucket, remote_dir)
            manifest = '\n'.join([
                f'WORKSPACE={remote}', 'BATCH_INDEX_OFFSET=0', 'CDAT=cdat.csv',
                'MAT=mat.fth', f'GROUP={group}', 'OUT_NAME=out.csv',
                f'MODEL=\'~group+n_genes{covs}\'', f'JOBS={jobs}',
            ])
            local_manifest = os.path.join(td, 'manifest.txt')
            with open(local_manifest, 'w') as m:
                m.write(manifest + '\n')
            remote_manifest = os.path.join(remote_dir, 'manifest.txt')
            print('Uploading manifest to s3...')
            s3.meta.client.upload_file(
                local_manifest, self.bucket, remote_manifest,
            )
        return remote_manifest

    def _mast_submit(
        self,
        manifest: str,
        block: bool = False,
        job_name: str = 'mast',
    ) -> str:
        batch = bt.client('batch')
        job_manifest = f's3://{os.path.join(self.bucket, manifest)}'
        job_id = None
        try:
            print(
                f'Submitting job {job_name} to the job queue {self.job_queue}'
            )
            submit_job_response = batch.submit_job(
                jobName=job_name, jobQueue=self.job_queue,
                jobDefinition=self.job_def,
                containerOverrides={'command': [job_manifest]}
            )
            job_id = submit_job_response['jobId']
            print(
                f'Submitted job {job_name} {job_id} to the job queue'
                f' {self.job_queue}'
            )
        except Exception as err:
            print(f'error: {str(err)}')
        return job_id

    def mast_compute(
        self,
        adata: AnnData,
        keys: Sequence[str],
        group: str,
        covs: str = '',
        block: bool = False,
        remote_dir: Optional[str] = None,
        jobs: int = 1,
    ) -> Tuple[str, str, str, Optional[DataFrame]]:
        content = None
        if remote_dir is None:
            remote_dir = os.path.join('mast', str(uuid.uuid4()))
            ready = []
        else:
            ready = ['mat', 'cdat']
        manifest = self._mast_prep(
            adata, remote_dir, keys, group, covs=covs, ready=ready, jobs=jobs,
        )
        job_name = f'mast-{"".join(filter(str.isalnum, group))}-{"".join(filter(str.isalnum, covs))}'
        job_id = self._mast_submit(
            manifest, block=block, job_name=job_name,
        )
        if block:
            status = BatchMAST._batch_job_status(job_id, wait=60)
            if status == 'SUCCEEDED':
                content = self._mast_results(remote_dir)
        return remote_dir, job_id, job_name, content

    def mast_collect(
        self,
        collection: Dict[str, Dict[str, str]],
        wait: int = 30,
    ) -> Generator[
        Tuple[str, str, Dict[str, str], Optional[DataFrame]],
        None,
        None,
    ]:
        while wait and len(collection) > 0:
            time.sleep(wait)
            for job_id in list(collection.keys()):
                status = BatchMAST._batch_job_status(job_id)
                if status == 'SUCCEEDED':
                    remote_dir = collection[job_id]['remote_dir']
                    content = self._mast_results(remote_dir)
                    yield job_id, status, collection.pop(job_id), content
                elif status == 'FAILED':
                    yield job_id, status, collection.pop(job_id), None

    @staticmethod
    def _batch_job_status(
        job_id: str,
        wait: int = 0,
        verbose: bool = False,
    ) -> str:
        batch = bt.client('batch')
        describe_jobs_response = batch.describe_jobs(jobs=[job_id])
        status = describe_jobs_response['jobs'][0]['status']
        if verbose:
            print(status)
        while wait:
            time.sleep(wait)
            describe_jobs_response = batch.describe_jobs(jobs=[job_id])
            new_status = describe_jobs_response['jobs'][0]['status']
            if new_status != status:
                status = new_status
                if verbose:
                    print(status)
            if status == 'SUCCEEDED' or status == 'FAILED':
                break
        return status

    def _mast_results(
        self,
        remote_dir: str,
    ) -> DataFrame:
        s3 = bt.resource('s3')
        with tempfile.TemporaryDirectory() as td:
            remote_out = os.path.join(remote_dir, 'out.csv')
            local_out = os.path.join(td, 'out.csv')
            s3.Bucket(self.bucket).download_file(remote_out, local_out)
            content = pd.read_csv(local_out, index_col=0)
        return content

    @staticmethod
    def mast_to_excel(
        de: Dict[str, DataFrame],
        fname: str,
        top: Optional[Dict[str, Dict[str, List[str]]]] = None,
    ):
        writer = pd.ExcelWriter(
            f'{fname}.xlsx', engine='xlsxwriter'
        )
        for s in de.keys():
            de[s].to_excel(writer, sheet_name=str(s))
        writer.save()
        if top is not None:
            writer = pd.ExcelWriter(
                f'{fname}.top.xlsx', engine='xlsxwriter'
            )
            for s in top.keys():
                pd.DataFrame.from_dict(
                    top[s], orient='index',
                ).T.fillna('').to_excel(
                    writer, sheet_name=str(s), index=False
                )
            writer.save()
