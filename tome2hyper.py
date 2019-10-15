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
    extract_table = TableDefinition(
        name=TableName("Reads", "Extract"),
        columns=[
            TableDefinition.Column(name='num_reads', type=SqlType.int(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name=sample_col, type=dtype(), nullability=NOT_NULLABLE),
            TableDefinition.Column(name=gene_col, type=dtype(), nullability=NOT_NULLABLE),
        ]
    )
    
    path_to_database = Path(hyper_file)        

    with HyperProcess(telemetry=Telemetry.SEND_USAGE_DATA_TO_TABLEAU) as hyper:

        # Replaces file with CreateMode.CREATE_AND_REPLACE if it already exists.
        with Connection(endpoint=hyper.endpoint,
                        database=path_to_database,
                        create_mode=CreateMode.CREATE_AND_REPLACE) as connection:

            connection.catalog.create_schema(schema=extract_table.table_name.schema_name)
            connection.catalog.create_table(table_definition=extract_table)

            # The rows to insert into the "Extract"."Extract" table.
            with Inserter(connection, extract_table) as inserter:
                for tome_file in tome_files:
                    tio = tomeio.TomeIO(tome_file)

                    all_samples = np.array(tio.sample_names)
                    all_genes = np.array(tio.gene_names)

                    if use_names: 
                        for si, ei, samples_index, gene_index, num_reads in tio.iter_gene_data(regions=regions):                    

                            sample_names = all_samples[samples_index]
                            gene_name = all_genes[gene_index]

                            for i in range(len(num_reads)):
                                inserter.add_row([int(num_reads[i]), sample_names[i], gene_name])
                    else:
                        for si, ei, samples_index, gene_index, num_reads in tio.iter_gene_data(regions=regions):
                            for i in range(len(num_reads)):
                                inserter.add_row([int(num_reads[i]), samples_index[i], gene_index])

                inserter.execute()

            # The table names in the "Extract" schema (the default schema).
            table_names = connection.catalog.get_table_names("Extract")
            print(f"Tables available in {path_to_database} are: {table_names}")

            # Number of rows in the "Extract"."Extract" table.
            # `execute_scalar_query` is for executing a query that returns exactly one row with one column.
            row_count = connection.execute_scalar_query(query=f"SELECT COUNT(*) FROM {extract_table.table_name}")
            print(f"The number of rows in table {extract_table.table_name} is {row_count}.")

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
    
    