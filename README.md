A Predictor of CD8+ T-cell epitope ImmunoGenicity given paired peptides and HLA-I alleles.

Briefly, PredIG score consists of a probability from 0 to 1, being 1 the
maximum likelihood for epitope immunogenicity. This score can be used to rank
candidates for prioritization approaches or for classification using adaptable
thresholds.

Please, refer to the Usage module for specifications of input format and data
requirements; refer to Method Specifics module for a detailed description on
feature space, model optimizations and user-choices.

# Tasks

## predig

Runs the complete pipeline, including the NOAH and NetCleave modules

### Input

PredIG requires as data points pairs of epitope and presenting HLA-I allele. 

For the file input, data points should be provided in a CSV with two columns named "Epitope" and "HLA_allele". CSV should be separated by commas. 

For the dialogue input, upload a list of epitope and HLA-I alleles in two separated columns entitled Epitope and HLA_allele, with commas as separators.

A job title can be assigned in the specific box. 

Please, load the example data for further clarifications.


###  Output

PredIG provides flexible tabular outputs including HTML, TSV and EXCEL.

The HTML output page can also be filtered and searched upon.

The results include the 14 features predicted within PredIG as well as an immunogenicity likelihood score termed PredIG score. 

Briefly, PredIG score consists of a probability from 0 to 1, being 1 the maximum likelihood for epitope immunogenicity. This score can be used to rank candidates for prioritization approaches or to classify them using adaptable thresholds. 

Please refer to Method Specifics section for a detailed description on feature space, model optimizations and user-choices.

