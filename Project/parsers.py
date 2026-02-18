from abc import ABC, abstractmethod
import pandas as pd


class FileParser(ABC):
    def __init__(self, file_path):
        self.file_path = file_path

    @abstractmethod
    def parse(self): #now every subclass of fileparser need to define a method called parse
        pass

class OBOParser(FileParser):

    def parse(self):
        rows = []
        current_term = None
        obsolete = False


        with open(self.file_path) as f:
            for line in f:
                line = line.strip()

                if line == "[Term]":
                    # save previous term
                    if current_term and not obsolete:
                        rows.append(current_term)

                    # start new term (ALWAYS)
                    current_term = {
                        "go_id": "",
                        "name": "",
                        "namespace": "",
                        "parents": [],
                        "definition":"",
                        "synonyms": []          #NEW  add this
                    }
                    obsolete = False #NEW reset


                elif current_term is not None:

                    if line.startswith("id:"):
                        current_term["go_id"] = line.split("id:")[1].strip()

                    elif line.startswith("name:"):
                        current_term["name"] = line.split("name:")[1].strip()

                    elif line.startswith("namespace:"):
                        current_term["namespace"] = line.split("namespace:")[1].strip()

                    elif line.startswith("is_a:"):
                        parent_id = line.split("is_a:")[1].split()[0]
                        current_term["parents"].append(parent_id)

                    elif line.startswith("is_obsolete: true"):
                        obsolete = True

                    elif line.startswith('def:'): # NEW Strips the OBO syntax around the definition and keeps only the content.
                        value = line.split('def:', 1)[1].strip()
                        if value.startswith('"'):
                            parts = value.split('"')
                            if len(parts) > 2:
                                current_term['definition'] = parts[1]
                            else:
                                current_term['definition'] = value
                        else:
                            current_term['definition'] = value

                    elif line.startswith("synonym:"): #NEW
                        # synonym: "TEXT" SCOPE [DBXREFS]
                        # we keep only the quoted TEXT
                        value = line.split("synonym:", 1)[1].strip()
                        if '"' in value:
                            # take the first quoted chunk
                            synonym_text = value.split('"')[1]
                            current_term["synonyms"].append(synonym_text)

        # save last term
        if current_term and not obsolete:
            rows.append(current_term)

        return pd.DataFrame(rows)



class GAFParser(FileParser):

    def parse(self):
        rows = []

        with open(self.file_path) as f:
            for line in f:
                if line.startswith("!"):
                    continue

                fields = line.strip().split("\t")

                rows.append({
                    "gene_id": fields[1],
                    "gene_name": fields[2],
                    'qualifier':fields[3],
                    "go_id": fields[4],
                    'aspect': fields [8],
                    "evidence": fields[6],
                    "molecule": fields[11]
                })

        return pd.DataFrame(rows)

