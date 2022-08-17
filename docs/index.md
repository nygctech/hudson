## Overview
hudson is a computational pipeline for spatial analysis of single cells in multiplexed fluorescence imaging. It is developed by Jagjit Singh, Kunal Pandit, and Rui Fu. The pipeline can installed by following the steps [here](https://github.com/nygctech/hudson)

Fluorescent Imaging, is the visualization of fluorescent dyes or proteins as labels for molecular processes or structures. It enables a wide range of experimental observations including the location and dynamics of gene expression, protein expression and molecular interactions in cells and tissues. fluorescence microscopy. 



IMAGE HERE


The pipeline has parallel capibility meaning it can be scaled over any cluster using the SLURM manager. It makes use of GPU nodes if avaiable for quick image processing and computation. The structure is defined in the next section.  

## Structure

The pipiline uses the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system to unite all its differnet components into one single workflow which can be replicated over any user machine or cluster. The different sections being: 

1) **Fix Lighting** : Explanation. 
2) **UNMIX**: Explanation. 
3) **Segmentation** : In the semgmenation step, the pipline performs instance segmentation on the image to generate a masks array of the same dimension as the input. The masks array is essentially the assignment of each pixel to a certain cell in the image. The users can select over which marker(LMN1b, GFAP, MBP, EVALV2,..) they want to segment on. Each cell is indexed as a label starting from 1 to the n'th cell that is segmented.
4) **Intensity Computation** : The next step is to compute the mean intensity per label using the masks array and the image. We need this to infer the cell types for the different labels that we got from the segmentation step. 
5) **Cell Type Indentification** 
6) **Data Object* 
7) **Spatial Neighborhood**:
8) 





### Markdown

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/nygctech/hudson/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.
