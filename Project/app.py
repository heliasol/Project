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
        qualifier=row["qualifier"],  
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

data = load_data()

template_data = {
    "total_terms": len(data["ontology_df"]),   #li prende da data, calc con miei metodi
    "leaf_perc": data["summary"]["leaf_percentage"],
    "avg_parents": data["summary"]["avg parents"],
    "avg_children": data["summary"]["avg children"],
    "exp_ratio": str((data['summary']['experimental vs computational'].get(True,0) / data['summary']['experimental vs computational'].sum() * 100)%),
    "mean_spec": data["specificity"]["mean_specificity"],
    "median_spec": data["specificity"]["median_specificity"],
    "genes_count": data["specificity"]["genes_analyzed"],
    "ns_bp": data["summary"]["namespace counts"].get("biological_process", 0),
    "ns_mf": data["summary"]["namespace counts"].get("molecular_function", 0),
    "ns_cc": data["summary"]["namespace counts"].get("cellular_component", 0),
    "ev_labels": data["summary"]["evidence counts"].index.tolist(),
    "ev_values": data["summary"]["evidence counts"].values.tolist(),
    "exp_count": data["summary"]["experimental vs computational"].get(True, 0),
    "comp_count": data["summary"]["experimental vs computational"].get(False, 0)
}

#routes
@app.route('/')
def home():
    return render_template('home.html')
    


@app.route("/term")
def term_page():
    go_id = request.args.get("go_id")

    term = None
    genes_for_term = []

    if go_id:
        term = terms.get_term(go_id)   # <-- fixed here
        if term:
            genes_for_term = annotations.get_by_term(go_id)

    return render_template(
        "term.html",
        go_id=go_id,
        term=term,
        genes_for_term=genes_for_term
    )
    


@app.route("/gene")
def gene_page():
    gene_name = request.args.get("gene_name")

    gene_annotations = []
    if gene_name:
        gene_annotations = annotations.get_by_gene_name(gene_name)
        print(f"Searching for gene: {gene_name}, found {len(gene_annotations)} annotations")

    return render_template(
        "gene.html",
        gene_name=gene_name,
        annotations=gene_annotations
    )

    

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
        
    return render_template ('Analyse_terms.html', result=result)

@app.route('/analyse_genes', methods= ['GET', 'POST'])
def analyse_genes():
    result= None
    if request.method == "POST":
        gene1 = request.form["gene1"]
        gene2 = request.form["gene2"]

        related = gene_analyser.genes_functionally_related(gene1, gene2)
        similarity_score= similarity_analyser.compare2genes(gene1,gene2)
        ancestor = gene_analyser.is_gene_ancestor(gene1, gene2)
        descendant = gene_analyser.is_gene_descendant(gene1, gene2)
        paths = gene_analyser.gene_paths(gene1,gene2)
        shortest_path = gene_analyser.shortest_gene_path(gene1, gene2)
        longest_path = gene_analyser.longest_gene_path(gene1, gene2)
        explanation = gene_analyser.explain_gene_relationship(gene1, gene2)

        result = {
            "gene1": gene1,
            "gene2": gene2,
            "related": related,
            'similarity_score': similarity_score,
            "ancestor": ancestor,
            "descendant": descendant,
            "path": paths,
            "shortest_path": shortest_path,
            'longest_path': longest_path,
            "explanation": explanation
        }

    return render_template("analyse_genes.html", result=result)

@app.route('/statistics')
def statistics():
    return render_template('statistics.html', **template_data)
    



