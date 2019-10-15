import pandas as pd
import numpy as np
import h5py
from memoized_property import memoized_property

class TomeIO:
    def __init__(self, file_name):
        self.f = h5py.File(file_name)                
            
    @memoized_property
    def samples(self):
        keys = list(self.f['sample_meta/anno'].keys())
        data_dict = {}
        for k in keys:
            dsv = self.f[f'sample_meta/anno/{k}'][()]
            if isinstance(dsv[0], bytes):
                data_dict[k] = [v.decode('utf-8') for v in dsv]
            else:
                data_dict[k] = dsv
                
        return pd.DataFrame(data=data_dict)            
        
    @memoized_property
    def gene_names(self):
        return [ g.decode("utf-8") for g in self.f["gene_names"][()] ]
    
    @memoized_property
    def sample_names(self):
        return [ s.decode("utf-8") for s in self.f["sample_names"][()] ]
    
    def iter_gene_data(self, regions, genes=None):
        all_genes = self.gene_names
        
        # h5py wants a boolean mask :/
        starts_index = np.ones(len(all_genes)+1, dtype=bool)
        starts_index[-1] = False
        ends_index = np.ones(len(all_genes)+1, dtype=bool)
        ends_index[0] = False
        
        if genes is not None:
            starts_index[:] = False
            ends_index[:] = False
            
            genes = set(all_genes) & set(genes)
            
            all_genes_index = np.array([ all_genes.index(g) for g in genes ])
            
            starts_index[all_genes_index] = True
            ends_index[all_genes_index+1] = True  
        else:
            all_genes_index = np.arange(len(all_genes))
                
        starts = self.f[f'data/{regions}/p'][starts_index]
        ends = self.f[f'data/{regions}/p'][ends_index]        
        
        start_i = 0        
        g_ct = 0
        for start, end, gidx in zip(starts, ends, all_genes_index):            
            g_ct += 1
            if g_ct % 2000 == 0:
                print(f'loaded genes {g_ct}/{len(all_genes_index)}') 
                
            num_reads = self.f[f'data/{regions}/x'][start:end]
            samples_index = self.f[f'data/{regions}/i'][start:end].astype(np.uint32)
            
            end_i = start_i + (end - start)
            
            yield start_i, end_i, samples_index, gidx, num_reads
            
            start_i = end_i 
        
    def read_gene_data_arrays(self, regions, genes=None):
                
        if genes is not None:
            num_rows = 0
            for start,end in zip(starts, ends):
                num_rows += end-start
        else:
            num_rows = self.f[f'data/{regions}/x'].shape[0]
            
        print(f'preallocating {num_rows} rows')
        
        out_data={
            'num_reads': np.empty(num_rows, dtype=np.uint32),
            'sample_index': np.empty(num_rows, dtype=np.uint32),
            'gene_index': np.empty(num_rows, dtype=np.uint32)
        }
        
        print("preallocated")
        
        for start_i, end_i, samples_index, gidx, num_reads in self.iter_gene_data(regions, genes):
            out_data['sample_index'][start_i:end_i] = samples_index
            out_data['gene_index'][start_i:end_i]= gidx
            out_data['num_reads'][start_i:end_i] = num_reads
            
        return out_data
    
    def read_gene_dataframe(self, regions, genes=None, use_names=False):
        data = self.read_gene_data_arrays(regions, genes)
        
        if use_names:
            print("converting index => name")
            data['sample_name'] = np.array(self.sample_names)[data['sample_index']]
            del data['sample_index']
        
            data['gene_name'] = np.array(self.gene_names)[data['gene_index']]
            del data['gene_index']
        
        print("creating dataframe")
        return pd.DataFrame.from_dict(data=data)                
    
    def read_gene_matrix(self, regions, genes=None, use_names=False):
        df = self.read_gene_dataframe(regions=regions, genes=genes, use_names=False)  
        
        df = df.pivot_table(index='sample_index', columns='gene_index', aggfunc=np.mean, fill_value=0.0)
        
        if use_names:            
            df.index = df.index.map(lambda idx: self.sample_names[idx])
            df.index.names = ['sample_name']
        
            df.columns = df.columns.map(lambda c: self.gene_names[c[1]])
            df.columns.names = ['gene_name']
        
        return df
    
    def __del__(self):
        self.f.close()