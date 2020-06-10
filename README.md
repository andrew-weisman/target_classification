# TARGET classification (using the GDC Data Portal)

## Directory setup

```bash
project_dir="/data/BIDS-HPC/private/projects/dmi2"
working_dir="/home/weismanal/notebook/2020-06-10/dmi"
mkdir "$project_dir" "$working_dir"
cd "$working_dir"
git clone git@github.com:andrew-weisman/target_classification.git "$project_dir/checkout"
mkdir "$project_dir/data"
```

Note: The effort using the data directly from the [TARGET data website](https://target-data.nci.nih.gov) (as opposed to the GDC Data Portal) is in the `target_data_website` branch of this repository.

## Workflow

Download the manifest for [all the gene expression quantification files in the TARGET program](https://portal.gdc.cancer.gov/repository?facetTab=files&files_size=100&files_sort=%5B%7B%22field%22%3A%22file_name%22%2C%22order%22%3A%22asc%22%7D%5D&filters=%7B%22op%22%3A%22and%22%2C%22content%22%3A%5B%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22cases.project.program.name%22%2C%22value%22%3A%5B%22TARGET%22%5D%7D%7D%2C%7B%22op%22%3A%22in%22%2C%22content%22%3A%7B%22field%22%3A%22files.data_type%22%2C%22value%22%3A%5B%22Gene%20Expression%20Quantification%22%5D%7D%7D%5D%7D&searchTableTab=files) (click on the blue "Manifest" button):

![all_gene_expression_files_in_target.png](images/all_gene_expression_files_in_target.png)

asdf
