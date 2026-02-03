## Isoform categories

SQANTI3 classifies each isoform by finding the best matching reference transcript, in the following order:

* **FSM (Full Splice Match)**: meaning the reference and query isoform have the same number of exons and each internal junction agree. The exact 5' start and 3' end can differ by any amount.

* **ISM (Incomplete Splice Match)**: the query isoform has fewer 5' exons than the reference, but each internal junction agree. The exact 5' start and 3' end can differ by any amount.

* **NIC (Novel In Catalog)**: the query isoform does not have a FSM or ISM match, but is using a combination of known donor/acceptor sites.

* **NNC (Novel Not in Catalog)**: the query isoform does not have a FSM or ISM match, and has at least one donor or acceptor site that is not annotated.

* **Antisense**: the query isoform does not have overlap a same-strand reference gene but is anti-sense to an annotated gene. 

* **Genic Intron**: the query isoform is completely contained within an annotated intron.

* **Genic Genomic**: the query isoform overlaps with introns and exons.

* **Intergenic**: the query isoform is in the intergenic region.

![sqanti_explain](https://raw.githubusercontent.com/FJPardoPalacios/public_figures/master/figuras_class_SQ3.png)

## Isoform subcategories

Some of the SQ3 categories are further divided into subcategories (specified in the `subcategory` field in SQANTI3 
classification output). These will be explained in the sections below.

### FSM subcategories

![FSM_subtype](https://raw.githubusercontent.com/FJPardoPalacios/public_figures/master/figure_fsm_subcat_SQ3.png)

### ISM subcategories

![ISM_subtype](https://raw.githubusercontent.com/FJPardoPalacios/public_figures/master/figure_ism_subcat_SQ3.png)

### Novel isoform subcategories (NIC and NNC)

Novel isoforms are subtyped based on whether they use a combination of known junctions (junctions are pairs of donor-acceptor sites), a combination of known splice sites (the individual donor and acceptor sites are known, but at least combination is novel), or at least one splice site (donor or acceptor) is novel.

![NIC_subtype](https://raw.githubusercontent.com/FJPardoPalacios/public_figures/master/figure_nic_nnc_subcat_SQ3.png)

