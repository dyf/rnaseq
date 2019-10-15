import argparse
from pathlib import Path
import tomeio
import numpy as np

from tableauhyperapi import HyperProcess, Telemetry, \
    Connection, CreateMode, \
    NOT_NULLABLE, NULLABLE, SqlType, TableDefinition, \
    Inserter, \
    escape_name, escape_string_literal, \
    TableName, \
    HyperException

def tome2hyper(tome_files, hyper_file, regions='exon', use_names=True):
    if use_names:
        sample_col = 'sample_name'
        gene_col = 'gene_name'
        dtype = SqlType.text()
    else:
        sample_col = 'sample_index'
        gene_col = 'gene_index'
        dtype = SqlType.int()
    
    dtype = SqlType.text if use_names else SqlType.int
    reads_table = TableDefinition(
        name=TableName("RNASeq", "Reads"),
        columns=[
            TableDefinition.Column(name='num_reads', type=SqlType.int(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name=sample_col, type=dtype(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name=gene_col, type=dtype(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name='tome_index', type=SqlType.int(), nullability=NOT_NULLABLE),
        ]
    )    

    tomes_table = TableDefinition(
        name=TableName("RNASeq", "Tomes"),
        columns=[
            TableDefinition.Column(name='tome_index', type=SqlType.int(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name='tome_name', type=SqlType.text(), nullability=NOT_NULLABLE)
        ]
    )

    samples_table = TableDefinition(
        name=TableName("RNASeq", "Samples"),
        columns=[
            TableDefinition.Column(name='sample_index', type=SqlType.int(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name='sample_name', type=SqlType.text(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name='tome_index', type=SqlType.int(), nullability=NOT_NULLABLE)
        ]
    )

    genes_table = TableDefinition(
        name=TableName("RNASeq", "Genes"),
        columns=[
            TableDefinition.Column(name='gene_index', type=SqlType.int(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name='gene_name', type=SqlType.text(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name='tome_index', type=SqlType.int(), nullability=NOT_NULLABLE)
        ]
    )

    path_to_database = Path(hyper_file)        

    with HyperProcess(telemetry=Telemetry.SEND_USAGE_DATA_TO_TABLEAU) as hyper:

        # Replaces file with CreateMode.CREATE_AND_REPLACE if it already exists.
        with Connection(endpoint=hyper.endpoint,
                        database=path_to_database,
                        create_mode=CreateMode.CREATE_AND_REPLACE) as connection:

            connection.catalog.create_schema(schema=reads_table.table_name.schema_name)
            connection.catalog.create_table(table_definition=reads_table)
            connection.catalog.create_table(table_definition=tomes_table)
            connection.catalog.create_table(table_definition=samples_table)
            connection.catalog.create_table(table_definition=genes_table)
            
            with Inserter(connection, samples_table) as inserter:
                for ti, tome_file in enumerate(tome_files):
                    tio = tomeio.TomeIO(tome_file)

                    for si, sample_name in enumerate(tio.sample_names):
                        inserter.add_row([si, sample_name, ti])
                inserter.execute()

            with Inserter(connection, genes_table) as inserter:
                for ti, tome_file in enumerate(tome_files):
                    tio = tomeio.TomeIO(tome_file)

                    for gi, gene_name in enumerate(tio.gene_names):
                        inserter.add_row([gi, gene_name, ti])
                inserter.execute()

            with Inserter(connection, tomes_table) as inserter:
                for ti, tome_file in enumerate(tome_files):
                    inserter.add_row([ti, tome_file])
                inserter.execute()

            print("wrote tomes")            
            with Inserter(connection, reads_table) as inserter:
                for tome_file in tome_files:
                    tio = tomeio.TomeIO(tome_file)

                    all_samples = np.array(tio.sample_names)
                    all_genes = np.array(tio.gene_names)

                    if use_names: 
                        for si, ei, samples_index, gene_index, num_reads in tio.iter_gene_data(regions=regions):                    

                            sample_names = all_samples[samples_index]
                            gene_name = all_genes[gene_index]

                            for i in range(len(num_reads)):
                                inserter.add_row([int(num_reads[i]), sample_names[i], gene_name, ti])
                    else:
                        for si, ei, samples_index, gene_index, num_reads in tio.iter_gene_data(regions=regions):
                            for i in range(len(num_reads)):
                                inserter.add_row([int(num_reads[i]), samples_index[i], gene_index, ti])

                inserter.execute()

            # The table names in the "Extract" schema (the default schema).
            table_names = connection.catalog.get_table_names("Extract")
            print(f"Tables available in {path_to_database} are: {table_names}")

            # Number of rows in the "Extract"."Extract" table.
            # `execute_scalar_query` is for executing a query that returns exactly one row with one column.
            row_count = connection.execute_scalar_query(query=f"SELECT COUNT(*) FROM {reads_table.table_name}")
            print(f"The number of rows in table {reads_table.table_name} is {row_count}.")

        print("The connection to the Hyper file has been closed.")
    print("The Hyper process has been shut down.")
    
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('hyper_file')
    parser.add_argument('tome_files', nargs='+')
    parser.add_argument('--use_names', action='store_true')
    parser.add_argument('--regions', default='exon')
    args = parser.parse_args()
    
    print(args)
    tome2hyper(args.tome_files, args.hyper_file, args.regions, args.use_names)
    
if __name__ == "__main__": main()
    
    