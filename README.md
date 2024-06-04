# Ciona Gene Model Converter

## About

This program converts gene identifiers between KY and KH gene models for _Ciona intestinalis_.  

Any issues can be reported to kpiekarz3@gatech.edu  
Image source (GUI): FABA https://www.bpni.bio.keio.ac.jp/chordate/faba/1.4/top.html  
Gene model source: http://ghost.zool.kyoto-u.ac.jp/default_ht.html  

The program originally coverts between KY21 and KH2013 gene models, but it also works with KH2012.  
There are two versions of the program: with and without the graphical user interface (GUI). GUI may be easier to run for an inexperienced user,
however, if you need to convert a large file, not just a few genes, it may be more convenient to use the version without GUI.

## How to get the code

To download the files, click on the green "Code" button above, select "Download ZIP", then unzip the folder.  

## How it works

Briefly, the equivalent genes are found by comparing their start positions (or end positions, depending on strandness). If there is an equivalent gene in the other gene model with the exact same start position, the program will return that gene. If there is no exact match, genes are sorted by the start position and the closest gene, either upstream or downstream, is returned. There is an adjustable parameter, called **distance cutoff** that sets the distance limit within which the program would search for an equivalent. 
#### Note on distance cutoff
Without a cutoff, the program will always find the closest gene, even if it's located hundreds of thousands of base pairs away, in which case it's most probably an unrelated gene, not a real equivalent. This is an arbitrary number, by default it's 500 bp, which was determined by trial and error to be a good starting point. In some cases, however, there are true equivalents with start positions even more than 50,000 bp apart, which can happen for example when a gene in KH gene model has just one exon, while its KY21 equivalent has for example 10 newly annotated exons (and it's therefore much longer), due to advances in sequencing/mapping, etc. So it's recommended to start with the default cutoff value, and only increase it for the queries that didn't produce a match during the initial search. However, keep in mind that increasing the cutoff can lead to errors, if another unrelated gene is located between the start positions of the true equivalents.

#### How to manually validate the results 
When in doubt, it is always a good idea to validate manually the results of the converter, at least for the key genes of interests.  
It can be done by going to http://ghost.zool.kyoto-u.ac.jp/default_ht.html.  
Choose the proper chromosome/contig from the dropdown list and paste the KY21 gene identifier in the search bar. When a list of transcriptional variants appears, choose any of them, and then check "KH2013 genes" on the leftside panel.  
<p align="center">
<img src="https://user-images.githubusercontent.com/117316002/205149747-cd032cc7-e3be-4484-bece-6010d01d7cf1.png" width="500">
</p>  

#### Possible outcomes, if no equivalent is found
* "No equivalent" - no equivalent found within the specified cutoff
* "No match" - when a contig where the query gene is located doesn't exist in the other gene model (a subtype of "No equivalent", used for troubleshooting)
* "No match in 2013" - this can happen when you try to convert a gene using a KH2012 identifier, when a gene is present in KH2012, but no longer exists in the updated, KH2013 gene model (it could have been merged with another gene, or is no longer considered a gene, etc.), as KH2012 identifier is subsequently compared to KH2013, before being converted to KY21

## How to run
* ### Version with GUI
There are several ways to run the program, some examples include:
1) Simply double-click on the executable file ("with_distance_cutoff_with_GUI_updated.exe") - no IDE or terminal use needed
2) Type "cmd" in the address bar in the file explorer (shown below), after the terminal opens, type "python with_distance_cutoff_with_GUI.py"  
<p align="center">
<img src="https://user-images.githubusercontent.com/117316002/205354321-5d181d24-83e1-41cc-8e5b-a6c41f63a356.png" width="270">
</p>  

3) Use any python IDE  

* ### Version without GUI
1) Open "with_distance_cutoff_noGUI.py" in a text editor and specify the input file (see "How to use" section below), then run it in a terminal as shown in point #2 above
2) Use directly in any python IDE

## How to use
* ### Version without GUI
Specify the input file in the "query" variable. For example:
```python
query = pd.read_csv("test_ky.csv")
```
The input should be a csv file with a column containing genes to be converted, with column name either "kh" or "ky".  
You can also modify the distance cutoff, if needed:
```python
cutoff = 2000
```
The output file ("converter_results.csv") will be generated automatically.
* ### Version with GUI
The input can be a single gene or a batch input (one gene per line).
Example input format: KY21.Chr1.1. (just gene identifier closed with a dot, without transcription variant annotations).
KY21 genes go into the left input field, the right input field is meant for KH genes. Paste your query in the appropriate input field and click "Get equivalent gene". Use the "Clear input" button to clear the input fields before submitting a new query.
The output is automatically copied to clipboard.  
<p align="center">
<img src="https://user-images.githubusercontent.com/117316002/205154763-eaa3bbf5-b3a6-4eca-94f3-de5a7b1f3c30.png" width="320">
</p>  

## Disclaimer

Since there is no 1:1 equivalency in gene positions in different gene models, sometimes the converter can erroneously identify a gene as an equivalent. Those mistakes can happen more often, if the distance cutoff is set to a large number, or if a gene is located in a problematic area with many overlapping genes. Always validate manually key genes of interest, before basing any decisions on the converter's output.

## How to cite

Katarzyna M Piekarz (2023). Ciona Gene Model Converter. Available from https://github.com/katarzynampiekarz/ciona_gene_model_converter
