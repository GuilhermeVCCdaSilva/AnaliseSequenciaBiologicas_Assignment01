#!/usr/bin/python
import toytree # Uma biblioteca para o desenho gráfico de árvores           
import toyplot.pdf


# Carrega um ficheiro newick com a arvore para a função do toytree
with open('bootstrap_MP12S.nhk','r') as file:
    bootstrap_MP12S = file.read()
tre = toytree.tree(bootstrap_MP12S)


# Uma lista com as diversas cores para os nomes das espécies da árvores
colorlist = ["darkcyan" if "Pfur" in t or "Lbra" in t or "Ltim" in t else "darkorange" for t in tre.get_tip_labels()]


# Define os estilos a utilizar
mystyle = {
    "edge_type": 'p',
    "edge_style": {
        "stroke": toytree.colors[2],
        "stroke": "darkcyan",
        "stroke": "darkorange",
        "stroke-width": 2.5,
    },
    "tip_labels_align": True,
    "tip_labels_colors" : colorlist,
    "tip_labels_style": {
        "font-size": "10px",
    },
    "node_labels" : tre.get_node_values("support",1,0),
    "node_sizes" : 7,
    "node_colors": toytree.colors[2],
}
# Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
canvas, axes, mark = tre.draw(height=700,**mystyle)


# Escreve a arvore em pdf
toyplot.pdf.render(canvas, "tree-plot_bootstrap_MP12S.pdf")


# Carrega um ficheiro newick com a arvore para a função do toytre
with open('bootstrap_NJCytb.nhk','r') as file:
    bootstrap_MP12S = file.read()
tre = toytree.tree(bootstrap_MP12S)

# Uma lista com as diversas cores para os nomes das espécies da árvores
colorlist = ["darkcyan" if "Pfur" in t or "Lbra" in t or "Ltim" in t else "indigo" for t in tre.get_tip_labels()]


# Define os estilos a utilizar
mystyle = {
    "edge_type": 'p',
    "edge_style": {
        "stroke": toytree.colors[2],
        "stroke-width": 2.5,
    },
    "tip_labels_align": True,
    "tip_labels_colors": colorlist,
    "tip_labels_style": {
        "font-size": "10px"
    },
    "node_labels" : tre.get_node_values("support",1,0),
    "node_sizes" : 7,
    "node_colors": toytree.colors[2],
}


# Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
canvas, axes, mark = tre.draw(height=900,**mystyle)


# Escreve a arvore em pdf
toyplot.pdf.render(canvas, "tree-plot_bootstrap_NJCytb.pdf")


# Carrega um ficheiro newick com a arvore para a função do toytree
with open('bootstrap_NJ12S.nhk','r') as file:
    bootstrap_MP12S = file.read()    
tre = toytree.tree(bootstrap_MP12S)

# Uma lista com as diversas cores para os nomes das espécies da árvores
colorlist = ["darkcyan" if "Pfur" in t or "Lbra" in t or "Ltim" in t else "darkorange" for t in tre.get_tip_labels()]


# Define os estilos a utilizar
mystyle = {
    "edge_type": 'p',
    "edge_style": {
        "stroke": toytree.colors[2],
        "stroke-width": 2.5,
        "stroke": "darkcyan",
        "stroke": "darkorange",
    },
    "tip_labels_align": True,
    "tip_labels_colors": colorlist,
    "tip_labels_style": {
        "font-size": "10px"
    },
    "node_labels" : tre.get_node_values("support",1),
    "node_sizes" : 7,
    "node_colors": toytree.colors[2],
}


# Guarda os valores que retornam do desenho da arvore em canvas, axes, mark
canvas, axes, mark = tre.draw(height=700,**mystyle)


# Escreve a arvore em pdf
toyplot.pdf.render(canvas, "tree-plot_bootstrap_NJ12S.pdf")


