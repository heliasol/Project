from flask import Flask, render_template, request
from parsers import OBOParser, GAFParser
from ontology import Term, TermCollection
from annotations import GeneAnnotation, AnnotationCollection
from hierarchy import OntologyHierarchy
from analysis import GeneAnalyser, GeneSimilarityAnalysis

app = Flask(__name__)

#parse the files
obo_df = OBOParser("go.obo").parse()
gaf_df = GAFParser("goa_human.gaf").parse()

#build ontology
terms = TermCollection()

for _, row in obo_df.iterrows():
    term = Term(
        go_id=row["go_id"],
        name=row["name"],
        namespace=row["namespace"],
        is_a=row["parents"],
        definition=row["definition"],
        synonyms=row["synonyms"]
    )
    terms.add_term(term)

terms.build_vertical_relationship()

#build annotations
annotations = AnnotationCollection()

for _, row in gaf_df.iterrows():
    ann = GeneAnnotation(
        gene_id=row["gene_id"],
        gene_name=row["gene_name"],
        go_id=row["go_id"],
        aspect=row["aspect"],
        evidence=row["evidence"],
        molecule=row["molecule"]
    )
    annotations.add_annotation(ann)

annotations.link_terms(terms)

#build hierarchy
hierarchy = OntologyHierarchy(terms)
hierarchy.build_tree()

#build analysers
gene_analyser = GeneAnalyser (annotations, terms, hierarchy)
similarity_analyser = GeneSimilarityAnalysis(obo_df, gaf_df)

#routes
@app.route("/")
def home():
    


@app.route("/term", methods=["GET", "POST"])
def term_page():
    


@app.route("/gene", methods=["GET", "POST"])
def gene_page():
    

@app.route("/analyse_terms", methods=["GET", "POST"])
def analyse_terms():
    result = None

    if request.method =='POST': 
        #runs the code only when the form is submitted
        #if the user is just opening the page normally
        #the code inside this block won't run
        go1 = request.form["go1"]
        go2 = request.form["go2"]

        related = hierarchy.is_related(go1,go2)
        ancestor = hierarchy.is_ancestor(go1, go2)
        descendant = hierarchy.is_descendant(go1,go2)
        msca = hierarchy.MSCA(go1,go2)
        paths = hierarchy.pedigree_paths(go1,go2)
        shortest_path=hierarchy.shortest_path(go1,go2)
        longest_path = hierarchy.longest_path(go1,go2)

        result = {
            'go1': go1,
            'go2' : go2,
            'related' : related,
            'ancestor' : ancestor,
            "descendant": descendant,
            "msca": msca,
            'paths': paths,
            "shortest_path": shortest_path,
            'longest_path': longest_path
            }
        
    return render_template ('compare_terms.html', result=result)

