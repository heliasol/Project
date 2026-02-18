from abc import ABC, abstractmethod
from parsers import *
from annotations import *
from ontology import *
from hierarchy import *
import numpy as np
import pandas as pd
from abc import ABC, abstractmethod
import matplotlib.pyplot as plt


class GeneAnalyser:
    def __init__(self,annotation_collection: AnnotationCollection,
                 term_collection: TermCollection,
                 hierarchy: OntologyHierarchy):

        self._annotations = annotation_collection
        self._ontology = term_collection
        self._hierarchy = hierarchy


    def _get_ann(self, gene: str) -> list[GeneAnnotation]: #it's private helper
        """
        Internal helper: return all annotations of a gene.
        """
        return self._annotations.get_by_gene_name(gene)
    
    def is_gene_ancestor(self, gene1: str, gene2: str) -> bool:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        for a1 in anns1:
            for a2 in anns2:
                if a1.term is not None and a2.term is not None:
                    if self._hierarchy.is_ancestor(a1.term.go_id,a2.term.go_id):
                        return True
        return False

    def is_gene_descendant(self, gene1: str, gene2: str) -> bool:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        
        for a1 in anns1:
            for a2 in anns2:
                if a1.term is not None and a2.term is not None:
                    if self._hierarchy.is_descendant(a1.term.go_id, a2.term.go_id):
                        return True
        return False
    
    def gene_specificity(self, gene: str) -> float | None:
        anns = self._get_ann(gene)
        depths = []
        
        for ann in anns:
            if ann.term is not None:
                ancestors = self._hierarchy._ontology.get_ancestors(ann.term.go_id)
                depths.append(len(ancestors))
        
        if not depths:
            return None

        return sum(depths) / len(depths)


    def genes_functionally_related (self, gene1: str, gene2: str) -> bool:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)

        terms1 = [ann.term.go_id for ann in anns1 if ann.term is not None]
        terms2 = [ann.term.go_id for ann in anns2 if ann.term is not None]


        for t1 in terms1:
            for t2 in terms2:
                if self._hierarchy.is_related(t1, t2):
                    return True
        return False


    def gene_paths(self, gene1: str, gene2: str):
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        
        paths = []

        for a1 in anns1:
            for a2 in anns2:
                if a1.term and a2.term:
                    subpaths = self._hierarchy.pedigree_paths(a1.term.go_id,a2.term.go_id)
                    paths.extend(subpaths)
        
        return paths
    
    def shortest_gene_path(self, gene1: str, gene2: str) -> list[str] | None:
        #highly related
        paths = self.gene_paths(gene1, gene2)
        return min(paths, key=len) if paths else None


    def longest_gene_path(self, gene1: str, gene2: str) -> list[str] | None:
       #more complex relationship
       paths = self.gene_paths(gene1, gene2)
       return max(paths, key=len) if paths else None
    
    def MSCA(self, gene1: str, gene2: str) -> str | None:
        anns1 = self._get_ann(gene1)
        anns2 = self._get_ann(gene2)
        
        for a1 in anns1:
            for a2 in anns2:
                if a1.term is None or a2.term is None:
                    continue
                
                msca = self._hierarchy.MSCA(a1.term.go_id, a2.term.go_id)

                if msca:
                    return msca   # return the first common ancestor found

        return None
