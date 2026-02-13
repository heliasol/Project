from ontology import *

class GeneAnnotation:
        branch_map = {
        "P": "Biological Process",
        "F": "Molecular Function",
        "C": "Cellular Component"} #NEW - fixed placement  # convert GAF aspect to branch

        def __init__(self, gene_id: str,        # DB Object ID (col 2)
        gene_name: str,      # DB Object Symbol (col 3)
        go_id: str,          # GO ID (col 5)
        qualifier: str | None = None,  # col 4 or relation in GAF 2.2
        aspect: str | None = None,     # col 9: P/F/C
        evidence: str | None = None,   # col 7
        molecule: str = "protein"  ):    # col 12: DB Object Type

            self._gene_id = gene_id
            self._gene_name = gene_name
            self._go_id = go_id
            self._qualifier = qualifier
            self._aspect = aspect
            self._evidence = evidence
            self._molecule = molecule

            self._term: Term | None = None       #NEW slot for Term






        @property
        def gene_name(self) -> str:
            """Human‑readable symbol from DB Object Symbol (column 3)."""
            return self._gene_name

        @property
        def go_id(self) -> str:
            """GO term identifier associated with this annotation."""
            return self._go_id

        @property
        def qualifier(self) -> str | None:
            """Qualifier / relation (e.g. 'NOT', 'contributes_to')."""
            return self._qualifier

        @property
        def aspect(self) -> str | None:
            """Aspect code P/F/C (BP/MF/CC)."""
            return self._aspect

        @property
        def branch(self) -> str:
            """Expanded aspect name (Biological Process, etc.)."""
            return self._branch

        @property
        def evidence(self) -> str | None:
            """Evidence code (IDA, IPI, IEA, etc.)."""
            return self._evidence

        @property
        def molecule(self) -> str:
            """DB Object Type, e.g. 'protein', 'gene', 'RNA'."""
            return self._molecule

        @property #
        def term(self) -> Term | None:
            """
            Linked Term object (if resolved); read‑only from outside.
            """
            return self._term

        def link_term(self, term_collection: "TermCollection") -> None:
            term = term_collection.get_term(self.go_id)
            if term != None:
                self._term = term


    #before go_id defined in init were just string
    #we want to replace it with Term object


        def __repr__(self) -> str:
            return f' gene name: {self.gene_name} \n go id: {self.go_id} \n branch: {self.branch} \n evidence: {self.evidence} \n molecule: {self.molecule} \n'


class AnnotationCollection: # holds all annotations and provides query functions

    def __init__(self) -> None:
        # internal mutable storage
        self._annotations: list[GeneAnnotation] = []

    @property
    def annotations(self) -> list[GeneAnnotation]:
        """
        Read-only style access: returns a shallow copy of the internal list.
        External code can iterate, filter, sort this copy without
        modifying the collection's internal state.
        """
        return self._annotations.copy()

    def __len__(self) -> int:
        return len(self._annotations)

    def __iter__(self):
        # Iterates over the internal list without exposing it directly.
        return iter(self._annotations)

    def __repr__(self) -> str:
        return f"<AnnotationCollection: {len(self._annotations)} annotations>"

    def add_annotation(self, annotation: GeneAnnotation) -> None:
        """Add a GeneAnnotation object to the collection."""
        self._annotations.append(annotation)

    def link_terms(self, term_collection: "TermCollection") -> None:
        """Resolve all GO IDs to Term objects."""
        for ann in self.annotations:
            ann.link_term(term_collection)

    def get_by_gene_id(self, gene_id: str) -> list[GeneAnnotation]:
        """
        Return all annotations for a specific gene_id (DB Object ID, col 2).
        """
        return [ann for ann in self._annotations if ann.gene_id == gene_id]

    def get_by_gene_name(self, gene_name: str) -> list[GeneAnnotation]:
        """
        Return all annotations for a specific gene symbol/name (col 3).
        """
        return [ann for ann in self._annotations if ann.gene_name == gene_name]

    def get_by_term(self, go_id: str) -> list[GeneAnnotation]:
        """Return all annotations for a specific GO term."""
        return [ann for ann in self._annotations if ann.go_id == go_id]

    def get_by_aspect(self, aspect: str) -> list[GeneAnnotation]:
        """
        Return all annotations with a specific aspect code (P, F, C).
        """
        return [ann for ann in self._annotations if ann.aspect == aspect]

    def get_by_evidence(self, evidence: str) -> list[GeneAnnotation]:
        """
        Return all annotations with a specific evidence code.
        """
        return [ann for ann in self._annotations if ann.evidence == evidence]

    def __repr__(self):
        return f"<AnnotationCollection: {len(self.annotations)} annotations>"