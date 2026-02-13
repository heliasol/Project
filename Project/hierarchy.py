from ontology import *
from annotations import *

class OntologyHierarchy:

    def __init__ (self, term_collection: TermCollection) -> None:
        self._ontology = term_collection
        self._hierarchy : dict[str, set[str]] = {}

    def build_tree(self) -> dict[str, set[str]]:
        for go_id in self._ontology.terms:
            self._hierarchy[go_id]= set()

        for term in self._ontology.terms.values():
            for parent in term.parents:
                self._hierarchy[parent.go_id].add(term.go_id)

        return self._hierarchy

    def is_descendant(self, child_id: str, parent_id: str) -> bool:
        descendants = self._ontology.get_descendants(parent_id)
        descendant_ids : set[str] = {t.go_id for t in descendants}
        return child_id in descendant_ids

    def is_ancestor(self, parent_id: str, child_id: str) -> bool:
        ancestors = self._ontology.get_ancestors(child_id)
        ancestor_ids : set[str] = {t.go_id for t in ancestors}
        return parent_id in ancestor_ids

    def is_related(self, go_id1: str , go_id2: str) -> bool:
        return (self.is_ancestor(go_id1,go_id2) or
                self.is_descendant(go_id1,go_id2))


    def pedigree_path(self, parent_id: str, child_id: str) -> list[str]:
      #I want to find a path from parent id to child id
      #output should look like parent id --> other ids if there --> child id (recursive)

        path = []

        if parent_id == child_id:
            return [parent_id]

        for child in self._hierarchy.get(parent_id, set()):
            p = self.pedigree_path(child, child_id)
            if p != []:
                return [parent_id] + p

        return []


    def pedigree_paths (self, parent_id: str, child_id: str) -> list[list[str]]:

        paths: list[list[str]] = []

        if parent_id == child_id:
            return [[parent_id]]

        for child in self._hierarchy.get(parent_id, set()):
            subpath = self.pedigree_paths(child, child_id)
            if subpath:
                for path in subpath:
                    paths.append([parent_id] + path)
        return paths

    def shortest_path(self, parent_id: str, child_id: str)  -> list[str] | None:
      paths = self.pedigree_paths(parent_id, child_id)
      if not paths:
            return None

      return min(paths, key=len) #key is a function that tells Python how to compare elements.


    def longest_path(self, parent_id: str, child_id: str)  -> list[str] | None:
      paths = self.pedigree_paths(parent_id, child_id)
      if not paths:
            return None
      return max(paths, key=len)


    def MSCA(self, go_id1: str, go_id2: str) -> str | None:
      #Most Specific Common Ancestor
      ancestors1 = self._ontology.get_ancestors(go_id1)
      ancestors2 = self._ontology.get_ancestors(go_id2)

      common = ancestors1.intersection(ancestors2)

      if not common:
        return None

      def depth(term) -> int:
        path = self.longest_path(term.go_id, go_id1)
        return len(path) if path is not None else 0

      return max(common, key=depth).go_id
    

    def __repr__(self):
        text =''
        for parent in self._hierarchy:

            if len(self._hierarchy[parent]) < 1:
                text += parent + '-->' + 'no children' + '\n'
            else:
                text += parent + '-->' + str(self._hierarchy[parent]) +'\n'

        return text