# JGV 1.0a2 package documentation

___
**JGV is an embed genomic viewer for Jupyter notebook written in python3**
___

# Installation

Ideally, before installation, create a clean python3 virtual environment to deploy the package, using virtualenvwrapper for example (see http://www.simononsoftware.com/virtualenv-tutorial-part-2/).

## Installation with pip from github

Install the package with pip3. All the required dependencies will be automatically installed.

```bash
pip3 install git+https://github.com/a-slide/JupyterGenoViewer.git
```

To update the package:

```bash
pip3 install git+https://github.com/a-slide/JupyterGenoViewer.git --upgrade
```

# Usage

The package is meant to be used in a jupyter notebook 4.0.0 +

## Notebook setup

Launch the notebook, navigate in the directory where you want to work and create a new python3 notebook

```bash
jupyter notebook
```

Import pylab (from matplotlib + numpy) and use %pylab magic command to enable plotting in the current Notebook.

```python
import pylab as pl
%pylab inline
```

    Populating the interactive namespace from numpy and matplotlib

Default pylab parameters can be defined at the beginning of the notebook as well (see http://matplotlib.org/users/customizing.html for more options)

```python
pl.rcParams['figure.figsize'] = 20,7
pl.rcParams['font.family'] = 'sans-serif'
pl.rcParams['font.sans-serif'] = ['DejaVu Sans']
pl.style.use('ggplot')
```

## Using JGV

JGV is first initialized with a reference genome. Then annotation and alignment files can be added. Finally, coverage and feature localization plots can be generated.

Each function has specific options that are comprehensively detailed in the testing notebook provided with the package or in html version on nbviewer: [Test_notebook](https://nbviewer.jupyter.org/github/a-slide/JupyterGenoViewer/blob/master/JGV/JGV_Test_Notebook.ipynb?flush_cache=true)

### Import package

```python
from JGV.JGV import JGV
```

One can also import the jprint and jhelp function from pycoQC to get a improve the default print and help function in jupyter

```python
from JGV.JGV import jhelp, jprint
```

A sample test file can be loaded from the package as well

```python
example_bam = JGV.example_bam()
example_fasta = JGV.example_fasta()
example_gtf = JGV.example_gtf()
example_gff3 = JGV.example_gff3()

jprint(example_bam)
jprint(example_fasta)
jprint(example_gtf)
jprint(example_gff3)
```

<p>/home/aleg/Programming/Python3/JupyterGenoViewer/JGV/data/yeast.bam</p>

<p>/home/aleg/Programming/Python3/JupyterGenoViewer/JGV/data/yeast.fa.gz</p>

<p>/home/aleg/Programming/Python3/JupyterGenoViewer/JGV/data/yeast.gtf.gz</p>

<p>/home/aleg/Programming/Python3/JupyterGenoViewer/JGV/data/yeast.gff3.gz</p>

### Initialize JGV with a reference genome

JGV starts by creating a Reference object from a fasta file

```python
j = JGV(fp=example_fasta, verbose=True)
```

<p><b>Add reference genome file</b></p>

<p>Parsing fasta file</p>

<p>&emsp;Found 17 reference sequences</p>

One can also give a list of chromosomes to select in the fasta file

```python
j = JGV(fp=example_fasta, verbose=True, ref_list=["I","II","III"])
```

<p><b>Add reference genome file</b></p>

<p>Parsing fasta file</p>

<p>&emsp;Found 17 reference sequences</p>

Finally, instead of a fasta file, one can provide a tab separated index file containing at least 2 columns with the refid(chromosome name) and the length of the sequence, such as a fasta index create by faidx or with the *output_index* option of JGV

```python
j = JGV(fp=example_fasta, verbose=True, output_index=True)
```

<p><b>Add reference genome file</b></p>

<p>Parsing fasta file</p>

<p>Write a fasta index file: /home/aleg/Programming/Python3/JupyterGenoViewer/JGV/data/yeast.tsv</p>

<p>&emsp;Found 17 reference sequences</p>

```python
index = "/home/aleg/Programming/Python3/JupyterGenoViewer/JGV/data/yeast.tsv"
j = JGV(index, verbose=True)
```

<p><b>Add reference genome file</b></p>

<p>Assume the file is a fasta index</p>

<p>&emsp;Found 17 reference sequences</p>

### Adding annotation files

Once initialized a JGV object can parse and save annotation files (gff3, gtf and bed).

```python
j.add_annotation(example_gtf, name="yeastMine")
```

Several annotation can be loaded. Warnings will be thrown if there are chromosomes found in the reference sequence have no feature in the annotation file

```python
j.add_annotation(example_gff3, name="Ensembl")
```

Information about the annotations can be obtained with annotation_summary

```python
j.annotation_summary()
```

<p><b>Counts per Annotation file</b></p>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Feature count</th>
      <th>Refid count</th>
      <th>Feature type count</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>yeastMine</th>
      <td>42071</td>
      <td>17</td>
      <td>6</td>
    </tr>
    <tr>
      <th>Ensembl</th>
      <td>28872</td>
      <td>17</td>
      <td>15</td>
    </tr>
  </tbody>
</table>

<p><b>Counts per Reference sequence</b></p>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>yeastMine</th>
      <th>Ensembl</th>
    </tr>
    <tr>
      <th>refid</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>I</th>
      <td>744</td>
      <td>509</td>
    </tr>
    <tr>
      <th>II</th>
      <td>2867</td>
      <td>1961</td>
    </tr>
    <tr>
      <th>III</th>
      <td>1174</td>
      <td>809</td>
    </tr>
    <tr>
      <th>IV</th>
      <td>5290</td>
      <td>3601</td>
    </tr>
    <tr>
      <th>IX</th>
      <td>1537</td>
      <td>1062</td>
    </tr>
    <tr>
      <th>Mito</th>
      <td>306</td>
      <td>257</td>
    </tr>
    <tr>
      <th>V</th>
      <td>2081</td>
      <td>1434</td>
    </tr>
    <tr>
      <th>VI</th>
      <td>916</td>
      <td>636</td>
    </tr>
    <tr>
      <th>VII</th>
      <td>3739</td>
      <td>2565</td>
    </tr>
    <tr>
      <th>VIII</th>
      <td>2027</td>
      <td>1391</td>
    </tr>
    <tr>
      <th>X</th>
      <td>2542</td>
      <td>1744</td>
    </tr>
    <tr>
      <th>XI</th>
      <td>2174</td>
      <td>1485</td>
    </tr>
    <tr>
      <th>XII</th>
      <td>3708</td>
      <td>2549</td>
    </tr>
    <tr>
      <th>XIII</th>
      <td>3236</td>
      <td>2226</td>
    </tr>
    <tr>
      <th>XIV</th>
      <td>2726</td>
      <td>1861</td>
    </tr>
    <tr>
      <th>XV</th>
      <td>3765</td>
      <td>2562</td>
    </tr>
    <tr>
      <th>XVI</th>
      <td>3239</td>
      <td>2220</td>
    </tr>
  </tbody>
</table>

<p><b>Counts per feature types</b></p>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>yeastMine</th>
      <th>Ensembl</th>
    </tr>
    <tr>
      <th>type</th>
      <th></th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>CDS</th>
      <td>7050.0</td>
      <td>7050.0</td>
    </tr>
    <tr>
      <th>chromosome</th>
      <td>NaN</td>
      <td>17.0</td>
    </tr>
    <tr>
      <th>exon</th>
      <td>7553.0</td>
      <td>7553.0</td>
    </tr>
    <tr>
      <th>gene</th>
      <td>7126.0</td>
      <td>6692.0</td>
    </tr>
    <tr>
      <th>mRNA</th>
      <td>NaN</td>
      <td>6692.0</td>
    </tr>
    <tr>
      <th>ncRNA_gene</th>
      <td>NaN</td>
      <td>15.0</td>
    </tr>
    <tr>
      <th>pseudogene</th>
      <td>NaN</td>
      <td>42.0</td>
    </tr>
    <tr>
      <th>rRNA</th>
      <td>NaN</td>
      <td>16.0</td>
    </tr>
    <tr>
      <th>rRNA_gene</th>
      <td>NaN</td>
      <td>16.0</td>
    </tr>
    <tr>
      <th>snRNA</th>
      <td>NaN</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>snRNA_gene</th>
      <td>NaN</td>
      <td>6.0</td>
    </tr>
    <tr>
      <th>snoRNA</th>
      <td>NaN</td>
      <td>77.0</td>
    </tr>
    <tr>
      <th>snoRNA_gene</th>
      <td>NaN</td>
      <td>77.0</td>
    </tr>
    <tr>
      <th>start_codon</th>
      <td>6700.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>stop_codon</th>
      <td>6516.0</td>
      <td>NaN</td>
    </tr>
    <tr>
      <th>tRNA_gene</th>
      <td>NaN</td>
      <td>299.0</td>
    </tr>
    <tr>
      <th>transcript</th>
      <td>7126.0</td>
      <td>314.0</td>
    </tr>
  </tbody>
</table>

### Adding alignment files

JGV objects can also parse and compute the coverage from alignment files (bam, sam and bed).

```python
j.add_alignment(example_bam, name="RNA-Seq")
```

Similar to annotation, JGV also has an alignment_summary function

```python
j.alignment_summary()
```

<p><b>Counts per Alignment file</b></p>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Refid count</th>
      <th>Base coverage</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>RNA-Seq</th>
      <td>17</td>
      <td>4051804</td>
    </tr>
  </tbody>
</table>

<p><b>Counts per Reference sequence</b></p>

<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>RNA-Seq</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>XII</th>
      <td>710501</td>
    </tr>
    <tr>
      <th>VII</th>
      <td>517868</td>
    </tr>
    <tr>
      <th>IV</th>
      <td>433392</td>
    </tr>
    <tr>
      <th>XV</th>
      <td>334195</td>
    </tr>
    <tr>
      <th>II</th>
      <td>271983</td>
    </tr>
    <tr>
      <th>XVI</th>
      <td>260855</td>
    </tr>
    <tr>
      <th>XI</th>
      <td>253186</td>
    </tr>
    <tr>
      <th>VIII</th>
      <td>236955</td>
    </tr>
    <tr>
      <th>X</th>
      <td>207338</td>
    </tr>
    <tr>
      <th>V</th>
      <td>203852</td>
    </tr>
    <tr>
      <th>XIII</th>
      <td>200794</td>
    </tr>
    <tr>
      <th>XIV</th>
      <td>130148</td>
    </tr>
    <tr>
      <th>III</th>
      <td>95877</td>
    </tr>
    <tr>
      <th>IX</th>
      <td>72737</td>
    </tr>
    <tr>
      <th>I</th>
      <td>72143</td>
    </tr>
    <tr>
      <th>VI</th>
      <td>49980</td>
    </tr>
    <tr>
      <th>Mito</th>
      <td>0</td>
    </tr>
  </tbody>
</table>

### Generate a plot of coverage per refid

Simple visualization to have a first idea of the sequencing coverage, with many customization options

```python
r = j.refid_coverage_plot()
```

![png](extra/extra/output_48_0.png)

```python
r = j.refid_coverage_plot(norm_depth=False, norm_len=False, log=True, color="dodgerblue", alpha=0.5)
```

![png](extra/extra/output_49_0.png)

### Plotting the coverage and annotation features of a specific window

interval_plot is undoubtedly the most useful function of the package. It has a large panel of option to customize the plots and will adapt automatically to plot all the annotation and alignment coverage over a defined genomic interval or an entire chromosome

```python
j.interval_plot("VI", feature_types=["gene", "transcript", "CDS"])
```

<p>Autodefine start position: 0</p>

<p>Autodefine end position: 270160</p>

<p>Estimated overlap offset: 675</p>

<p><b>Extract alignment data</b></p>

<p>Compute coverage from the windows: VI:0-270160</p>

<p>&emsp;Define size of each bin: 540.32</p>

<p>&emsp;Compute coverage...</p>

<p><b>Extract annotation data</b></p>

<p>&emsp;Alignment track name: RNA-Seq</p>

<p>&emsp;Alignment track name: yeastMine</p>

<p>&emsp;Alignment track name: Ensembl</p>

![png](extra/extra/output_52_11.png)

```python
j.interval_plot("VI", start=220000, end=225000)
```

<p>Estimated overlap offset: 12</p>

<p><b>Extract alignment data</b></p>

<p>Compute coverage from the windows: VI:220000-225000</p>

<p>&emsp;Define size of each bin: 10.0</p>

<p>&emsp;Compute coverage...</p>

<p><b>Extract annotation data</b></p>

<p>&emsp;Alignment track name: RNA-Seq</p>

<p>&emsp;Alignment track name: yeastMine</p>

<p>&emsp;Alignment track name: Ensembl</p>

![png](extra/extra/output_53_9.png)

# Note to developers

You are welcome to contribute by requesting additional functionalities, reporting bugs or by forking, and submitting a pull request.

Thank you

# Authors and Contact

Adrien Leger - 2017

Enright's group, EMBL EBI

* <aleg@ebi.ac.uk>
* [Github](https://github.com/a-slide)
