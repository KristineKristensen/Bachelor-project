Scripts for inferring a species tree for Ceroxyloideae with paralogs

1. Programs used
   
    - HybPiper 2.1.2 (Johnson et al., 2016)
    - BWA-MEM 0.7.17 (Li, 2013) 
    - Picard Markdeduplicate 2.27.5 (Broad institute)
    - Samtools 1.3.1. (Danecek et al., 2021)
    - MAFFT version 7.520 (Katoh and Standley, 2013)
    - TrimAl version 1.4.rev15 (Capella-Gutiérrez et al., 2009)
    - AMAS 1.0 (Borowiec, 2016)
    - R-4.3.0 0 (R Core Team, 2023) 
    - CIAlign 1.0.18 (Tumescheit et al., 2022)
    - TAPER Version 1.0.0 (Zhang et al., 2021) 
    - IQtree multicore version 2.2.2.3 (Minh et al., 2020) with ModelFinder Plus (Kalyaanamoorthy et al., 2017) 
    - ASTRAL-III version 5.7.8 (Zhang et al., 2018) 
    - Phytools 1.5-1 (Revell, 2012)
    
  
 
 
 


2. Analysis

All analyses were carried out on genomeDK following Workflow.py. All scripts was provided by Paola de Lima Ferreira.
 
1. HypbPiper assemble
2. HybPiper stats
3. HybPiper heatmap
4. HybPiper paralog
5. Coverage, at this step is the script coverage.py is used (a version of https://github.com/Sarah-Elna/BSc/blob/2d644b856ee34ece75355dcf5489bc846bee2453/Python_Scripts/coverage_eddit.py)
7. Retrieve, at this step the sample2genes.py was used (from https://github.com/pebgroup/Dypsidinae_species_tree/blob/16c8082baad04b8d5e7efc2052aadc132660bae0/samples2genes.py).
7. Mafft
8. Trimal
9. AMAS, for this step the scripts amas_raw.sh and amas_gt.sh was used
11. Optimal, here the Optrimal.R script was run (a version of the script from Shee Z.Q., Frodin D.G., Cámara-Leret R. and Pokorny L. 2020. Reconstructing the Complex Evolutionary History of the Papuasian Schefflera Radiation Through Herbariomics. Frontiers in Plant Science. 11. https://www.frontiersin.org/articles/10.3389/fpls.2020.00258)
12. Cialign
13. Taper
14. IQtree
15. Astral
	


