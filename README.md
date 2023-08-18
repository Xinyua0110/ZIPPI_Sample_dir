# ZIPPI_log
1. using the Sample_dir provided by Prof.Haiqing Zhao to **caculate the 6 metrics** of selected 100 bacteria PDB dimers, (50 from bacteria PDB heterodimers, and 50 from bacteria PDB homodimers),[coding](https://github.com/Xinyua0110/ZIPPI_Sample_dir/tree/6-metrics)
   * mutual information (MI), 
   * conservation (Con) and 
   * direct coupling (DCA). 
   * average product correction (APC) of MI
   * average product correction (APC) of Con
   * average product correction (APC) of DCA

2. Test the caculation on CAPRI benchmark decoys-[score_set](https://scoreset.org/index.php?browse)

   * Caculate interface residues : saved in xxx.ifc_seq file

     [Example](https://github.com/Xinyua0110/ZIPPI_Sample_dir/blob/score_set/get_interface_contacts.py)

     ‘The iRMSD is the Cα RMSD of interface residues with respect to the native structure. <u>Interface residues in a complex are defined as all the residues within 10.0 Å from any residues of the other subunit</u>.’——Wang, X., Flannery, S. T., & Kihara, D. (2021). Protein Docking Model Evaluation by Graph Neural Networks. *Frontiers in Molecular Biosciences*, *8*, 647915. https://doi.org/10.3389/fmolb.2021.647915

   * The msa file-paired based on the shared common species.:save in xxx.msa

     To avoid biased sequence sampling due to over-studied model species, we carried out homolog sequence search on 5,090 representative proteomes that were carefully curated and selected in **EggNog 5.0.**

     `/home/zjlab/score_set/e5.proteomes.faa:`http://eggnog5.embl.de/download/eggnog_5.0/e5.proteomes.faa

     `/home/zjlab/score_set/Q03774.fasta :`https://www.ebi.ac.uk/pdbe/pdbe-kb/proteins/Q03774

     ```
     jackhmmer -N 5 --incE 0.001 --tblout output.tblout /home/zjlab/score_set/Q03774.fasta /home/zjlab/score_set/e5.proteomes.faa > output.sto
     ```

     

   * Use the two above to caculate 6 metrics

   * compare the caculated 6 metrics to the ZIPPI paper -SI Table S2. Detailed AUROC performances of ZIPPI metrics for each CAPRI Target
