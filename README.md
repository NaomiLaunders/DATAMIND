### DATAMIND: Roadbuilder 3: Investigating severe mental health and physical illness across the four nations

Datamind is a collaborative project aiming to enrich mental health data science in the UK. I am working on RoadBuilder 3, which seeks to catalogue the routine data resources available in England, Wales, Scotland and Northern Ireland for studying severe mental illness and physical health. 
[Read more about Datamind here!](https://www.hdruk.ac.uk/helping-with-health-data/health-data-research-hubs/datamind/)

We tested the feasibility of studying cardiovascular disease risk factors in people with severe mental illness across a range of data sources. 

Here are our reports: https://datamind.org.uk/data/harmonised-data/linking-mental-and-physical-health/

And two publications:

[Prevalence and patient characteristics associated with cardiovascular disease risk factor screening in UK primary care for people with severe mental illness: an electronic healthcare record study](https://pubmed.ncbi.nlm.nih.gov/39819835/)

[Characteristics of people with severe mental illness excluded from incentivised physical health checks in the UK: electronic healthcare record study](https://pubmed.ncbi.nlm.nih.gov/40377164/)

## Understanding the files
Sample scripts have been uploaded to this GitHub page. As a first use of GitHub and of looking at some of these measures in CPRD, I will be optimising some of this code for future studies!

Scripts are organised into chunks: Code lists (scripts beginning 1), observations (scripts beginning 2), cohort creation (scripts beginning 3) and analysis (scripts beginning 4).

Only a subset of scripts have been provided.

Code lists are available on the [HDR-UK phenotype library](https://phenotypes.healthdatagateway.org/)  

Code lists were generated in CPRD. They are specific to CPRD and to the exact database builds. When re-using codelists it's advisable to search for new terms in the system being used.

CPRD code lists use medcodes. These long numeric identifiers must be kept as character vectors in R, Stata and Excel to avoid truncation.

This project uses pseudonymised patient records and so no data is available for uploading.
