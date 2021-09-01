##EFO and MONDO mapping

The Experimental Factor Ontology (EFO) is a resource that provides systematic description of experimental variables including cancer groups available in EBI databases and some other external projects. The scope of EFO is to support the annotation, analysis and visualization of data handled by many groups at the EBI and as the core ontology for OpenTargets.org


The Mondo Disease Ontology (Mondo) is another independent resource aiming to harmonize disease definitions. 

The purpose of maintaining `efo-mondo-map.tsv` is to map each cancer group found in the `histologies.tsv` to its appropriate EFO and MONDO codes.

###Format

The EFO codes may have prefix “EFO_”, “MONDO_” or “Orphanet_”, while the MONDO codes have a prefix “MONDO_”; each followed by seven digits.

As more cancer groups or subtypes are added to existing histologies dataset and/or in case of any ambiguity, these codes will be revised.
Currently, all non-NA cancer groups from `histologies.tsv` except 'Epilepsy', 'Dysplasia/Gliosis;Glial-neuronal tumor NOS', 'Arteriovenous malformation', and 'Reactive connective tissue' are included `efo-mondo-map.tsv`.
