# HappySorter pipeline

As preparation, you should have datastores for PE and long reads. Also you should know the unique 31-mer coverage peak for the PE, but you allways do a 31-mer spectra on your PE datasets anyway, don't you?

This pipeline goes in four simple steps:

* **01-dbg_strider.py** creates a dbg and extend its haplotype-specific unitigs with short and long reads.

* **02-anchors_rtg.py** Create anchors and a ReadThreadsGraph.

* **03-make_orders.py** Uses the TotalSorter/HappySorters to create and merge the orders. Includes post-order scaffolding.

* **04-make_outputs.py** Creates output files (assembly consensus, etc).