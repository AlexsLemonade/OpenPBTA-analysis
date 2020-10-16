# Molecular Subtyping of Neurocytoma

**Module authors:** Krutika Gaonkar ([@kgaonkar6](https://github.com/kgaonkar6i))

In this analysis we subtype Neurocytoma samples according to the primary_site and pathology_diagnosis values. If primary_site == "Ventricles" then subtype is CNC (central neurocytoma) if primary_site != "Ventricles" then subtype is EVN (extraventricular neurocytoma)

### Inclusion/exclusion criteria

Samples are _included_ for subtyping if we detect the following strings in the `pathology_diagnosis` field of `pbta-histologies.tsv`:

```
Neurocytoma
```

`01-neurocytoma-subtyping.Rmd` reads in the histology file and adds molecular_subtype for neurocytoma samples as `CNC` if primary_site == "Ventricles" and `EVN` if primary_site != "Ventricles"  

##o Run script

```sh
bash run_subtyping.sh
```

