


<!-- PROJECT LOGO -->
<br />
<p align="center">
  <a href=>
  <h3 align="center"> Balanced SET levels favor the correct enhancer repertoire during cell fate acquisition
</h3>
  </a>
  <p align="center">
    Mattia Zaghi<sup>1 *</sup>, Federica Banfi<sup>1</sup>, Luca Massimino<sup>1</sup>, Monica Volpin<sup>1</sup>, Edoardo Bellini<sup>1</sup>, Simone Brusco<sup>1</sup>, Ivan Merelli<sup>1</sup>, Cristiana Barone<sup>1</sup>, Cristina Sironi<sup>1</sup>, Michela Bruni<sup>1</sup>, Linda Bossini<sup>1</sup>, Luigi Lamparelli<sup>1</sup>, Laura Pintado<sup>1</sup>, Deborah D'Aliberti<sup>1</sup>, Luca Mologni<sup>1</sup>, Gaia Colasante<sup>1</sup>, Federica Ungaro<sup>1</sup>, Francesco Ferrari<sup>1</sup>, Jean-Michel Cioni<sup>1</sup>, Emanuele Azzoni<sup>1</sup>, Rocco Piazza<sup>1</sup>, Eugenio Montini<sup>1</sup>, Vania Broccoli<sup>1</sup> and Alessandro Sessa<sup>1</sup>
  </p>
</p>
<p align="center">
<img src="./image/UMAP.jpeg" width="300" height="300" />
<img src="./image/UMAP_ATAC.jpeg" width="300" height="300" />
</p>
<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary>Table of Contents</summary>
  <ol>
  <li>
      <a href="#About the paper">About the paper</a>
  </li>
    <li>
      <a href="#pipelines">Pipelines</a>
    </li>
    <li>
      <a href="#figures">Figures</a>
    </li>
  </ol>
</details>



<!-- About the paper -->
## About the paper

In this work we described the effect of SET protein accumulation on chromatin rewiring in vitro, using Schinzel-Giedion syndrome patients IPSCs, and in vivo using Mouse and Zebrafish model. We used a multiomic approach combining ATAC-seq, ChIP-seq, Hi-C, RNA-seq and scMultiome (ATAC+RNA) to address the question on a genomic standpoint. In this page all the codes and pipelines used to analyze all data and produce the manuscript figures are deposited.


<!-- Pipelines -->
## Pipelines
In this folder the snakemake pipelines, bash and R script used to analyzed all the epigenomic/trascriptomic data (ATAC, ChIP-seq, Hi-C, RNA-seq, scMultiome) are deposited. Within the folder are present:

* `.sk` files, that contain all rules sequence.

* `.yaml` file, that must be edit according to organism and statistics and user needs.

* `.R` files, custom R script use to analyze a specific dataset.

<!-- Figures -->
## Figures
In Figures folder are present all the costum snakemake and `.R` script used to generate plots and figures. Each figure folder contain the code used to generate it and its related supplementary figure.

<!-- CONTACT -->
## Contacts

Mattia Zaghi - <zaghi.mattia@hsr.it>

Edoardo Bellini - <bellini.edoardo@hsr.it>

Alessandro Sessa - <sessa.alessandro@hsr.it>

<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/othneildrew/Best-README-Template.svg?style=for-the-badge
[contributors-url]: https://github.com/othneildrew/Best-README-Template/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/othneildrew/Best-README-Template.svg?style=for-the-badge
[forks-url]: https://github.com/othneildrew/Best-README-Template/network/members
[stars-shield]: https://img.shields.io/github/stars/othneildrew/Best-README-Template.svg?style=for-the-badge
[stars-url]: https://github.com/othneildrew/Best-README-Template/stargazers
[issues-shield]: https://img.shields.io/github/issues/othneildrew/Best-README-Template.svg?style=for-the-badge
[issues-url]: https://github.com/othneildrew/Best-README-Template/issues
[license-shield]: https://img.shields.io/github/license/othneildrew/Best-README-Template.svg?style=for-the-badge
[license-url]: https://github.com/othneildrew/Best-README-Template/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/othneildrew
