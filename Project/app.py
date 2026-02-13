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
    return render_template("index.html")


@app.route("/term", methods=["GET", "POST"])
def term_page():
    result = None

    if request.method == "POST":
        go_id = request.form["go_id"]
        term = terms.get_term(go_id)

        if term:
            result = {
                "id": term.go_id,
                "name": term.name,
                "namespace": term.namespace,
                "definition": term.definition,
                "parents": [p.go_id for p in term.parents],
                "children": [c.go_id for c in term.children]
            }

    return render_template("term.html", result=result)


@app.route("/gene", methods=["GET", "POST"])
def gene_page():
    summaries = None
    specificity = None
    gene = None

    if request.method == "POST":
        gene = request.form["gene_name"]

        anns = annotations.get_by_gene_name(gene)

        summaries = [
            gene_analyser.summary(a)
            for a in anns
            if gene_analyser.summary(a)
        ]

        specificity = gene_analyser.gene_specificity(gene)

    return render_template("gene.html",
                           gene=gene,
                           summaries=summaries,
                           specificity=specificity)


@app.route("/compare", methods=["GET", "POST"])
def compare_page():
    related = None
    similarity_score = None
    explanation = None

    if request.method == "POST":
        gene1 = request.form["gene1"]
        gene2 = request.form["gene2"]

        related = gene_analyser.genes_functionally_related(gene1, gene2)
        explanation = gene_analyser.explain_gene_relationship(gene1, gene2)

        similarity_score = similarity_analyser.compare2genes(gene1, gene2)

    return render_template("compare.html",
                           related=related,
                           explanation=explanation,
                           similarity=similarity_score)


if __name__ == "__main__":
    app.run(debug=True)
