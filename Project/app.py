from flask import Flask, render_template, request
from parsers import OBOParser, GAFParser
from ontology import Term, TermCollection
from annotations import GeneAnnotation, AnnotationCollection
from hierarchy import OntologyHierarchy
from analysis import *

# initialize app
app = Flask(__name__)

def load_data():
    # parse the files
    print('parsing start')
    obo_df = OBOParser("gene ontology.txt").parse()
    gaf_df = GAFParser("gaf.txt").parse()
    print('parsing done')

    # build ontology
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

    print('term relationship is being made')
    terms.build_vertical_relationship()
    print('term relationship made')

    # build annotations
    print('annotation starts')
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

    print('annotation made')

    print('linking on process')
    annotations.link_terms(terms)
    print('links made')

    # build hierarchy
    print('building hierarchy')
    hierarchy = OntologyHierarchy(terms)
    hierarchy.build_tree()
    print('hierarchy tree done')

    # build analysers
    print('start analysing')
    gene_analyser = GeneAnalyser(annotations, terms, hierarchy)
    print('finish analysing')

    #stat
    print('start summary')
    summary= SummaryStatistics(obo_df, gaf_df, terms)
    print('finish summary')

    #similarity analysis
    print('similaity start')
    similarity_analyser=GeneSimilarityAnalysis(obo_df,gaf_df)
    print('similarity finish')

    # return structured data
    return {
        "ontology_df": obo_df,
        "annotation_df": gaf_df,
        "term_collection": terms,
        'annotations': annotations,
        'hierarchy': hierarchy,
        "gene_analyser": gene_analyser,
        "summary": summary,
        'similarity_analyser':similarity_analyser}
    


data = load_data()
hierarchy=data['hierarchy']
gene_analyser=data['gene_analyser']
terms=data['term_collection']
annotations=data['annotations']
similarity_analyser = data['similarity_analyser']
ontology_df= data['ontology_df']
annotation_df= data ['annotation_df']



    
#routes
@app.route('/')
def home():
    return render_template('index.html')

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


@app.route("/term")
def term_page():
    go_id = request.args.get("go_id")

    term = None
    genes_for_term = []

    if go_id:
        term = terms.get_term(go_id)  
        if term:
            genes_for_term = annotations.get_by_term(go_id)

    return render_template(
        "term.html",
        go_id=go_id,
        term=term,
        genes_for_term=genes_for_term
    )
    

@app.route("/analyse_terms", methods=["GET", "POST"])
def analyse_terms():
    result = None

    if request.method =='POST': 
        go1 = request.form["go1"]
        go2 = request.form["go2"]

        result = {
            'go1': go1,
            'go2' : go2,
            'related': hierarchy.is_related(go1,go2),
            'ancestor' : hierarchy.is_ancestor(go1, go2) ,
            "descendant": hierarchy.is_descendant(go1,go2),
            "msca": hierarchy.MSCA(go1,go2),
            'paths': hierarchy.pedigree_paths(go1,go2),
            'altpaths':hierarchy.pedigree_paths(go2,go1),
            "shortest_path": hierarchy.shortest_path(go1,go2),
            'altshortest_path': hierarchy.shortest_path(go2,go1),
            'longest_path':hierarchy.longest_path(go1,go2),
            'altlongest_path': hierarchy.longest_path(go2,go1)
            }

 
    return render_template ('analyse_terms.html',
                            result=result ,
                            go1=result['go1'] if result else None,
                            go2=result['go2'] if result else None)


@app.route('/analyse_genes', methods= ['GET', 'POST'])
def analyse_genes():
    result= None
    if request.method == "POST":
        gene1 = request.form["gene1"]
        gene2 = request.form["gene2"]

        result = {
            "gene1": gene1,
            "gene2": gene2,
            'similarity_score': similarity_analyser.compare2genes(gene1,gene2),
            "ancestor": gene_analyser.is_gene_ancestor(gene1, gene2),
            "descendant": gene_analyser.is_gene_descendant(gene1, gene2),
            "paths": gene_analyser.gene_paths(gene1,gene2),
            'altpaths': gene_analyser.gene_paths(gene2,gene1),
            "shortest_path": gene_analyser.shortest_gene_path(gene1, gene2),
            'longest_path': gene_analyser.longest_gene_path(gene1, gene2),
            'altshortest_path': gene_analyser.shortest_gene_path(gene2,gene1),
            'altlongest_path':gene_analyser.longest_gene_path(gene2,gene1),
            'msca': gene_analyser.MSCA(gene1,gene2)
        }

    return render_template("analyse_genes.html", result=result, gene1=result['gene1'] if result else None, gene2=result['gene2'] if result else None)

@app.route('/stats')
def stats(): 
    return 'Hello for now'



if __name__ == '__main__':
    app.run(debug=True, use_reloader=False)
