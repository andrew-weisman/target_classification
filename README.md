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

asdf
