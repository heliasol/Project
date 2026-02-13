class Term: #Represents a single GO term.(nodes)  #property
    def __init__(self, go_id: str, name: str, namespace: str,
                 is_a: list[str] | None, definition: str,
                 synonyms: list[str] | None = None):
        self.__go_id = go_id
        self.__name = name
        self.__namespace = namespace
        self.__is_a = is_a or []  # raw parent IDs from parsing
        self.__definition = definition
        self.__synonyms = synonyms or []


        self.__parents = set()       # Term objects
        self.__children = set()      # Term objects

#we use copy for sets, lists and dicts (mutable) and we don't need to use it for str and int(immutable)

    @property
    def go_id(self):
        return self.__go_id

    @property
    def name(self):
        return self.__name

    @property
    def namespace(self):
        return self.__namespace

    @property
    def definition(self):
        return self.__definition

    @property
    def synonyms(self):
        return self.__synonyms

    @property
    def parents(self):
        return self.__parents

    @property
    def children(self):
        return self.__children

    @property
    def is_a(self):
        return self.__is_a


    def add_parent(self, parent: 'Term'):
        self.__parents.add(parent)
        parent.__children.add(self)

    def __repr__(self):
        return f" *GOTerm \n go id: {self.__go_id} \n name: {self.__name} \n namespace: {self.__namespace} \n parents: {self.__parents} \n children: {self.__children} \n"
#Now we have full collection of terms => we can start interconnecting them


class TermCollection: #the entire graph
    def __init__(self) -> None:
        self.__terms: dict[str, Term] = {}

    @property #Lia added this
    def terms(self) -> dict[str, Term]:
      return self.__terms

    def add_term(self, term: Term) -> None:
        self.__terms[term.go_id] = term


    def build_vertical_relationship(self): # it creates relationships from the is_a thing
            for term in self.__terms.values():
                for parent_id in term.is_a:
                    parent = self.__terms.get(parent_id)
                    if parent != None:
                        term.add_parent(parent)


    def get_term(self, go_id: str) -> Term | None:
            return self.__terms.get(go_id)


    def get_parents(self, go_id: str) -> set[Term]:
        term = self.get_term(go_id)
        return term.parents if term else set()

    def get_children(self, go_id: str) -> set[Term]:
        term = self.get_term(go_id)
        return term.children if term else set()

    def get_ancestors(self, go_id: str) -> set[Term]:
            """Return all ancestor terms of a given GO term (recursively)."""
            term = self.get_term(go_id)
            if not term:
                return set()

            ancestors: set[Term] = set()

            def explore(term: Term) -> None:
                for parent in term.parents:
                    if parent not in ancestors:
                        ancestors.add(parent)
                        explore(parent)

            explore(term)
            return ancestors

    def get_descendants(self, go_id: str) -> set[Term]:
        """Return all descendant terms of a given GO term (recursively)."""
        term = self.get_term(go_id)
        if not term:
            return set()

        descendants: set[Term] = set()

        def explore(term: Term) -> None:
            for child in term.children:
                if child not in descendants:
                    descendants.add(child)
                    explore(child)

        explore(term)
        return descendants
