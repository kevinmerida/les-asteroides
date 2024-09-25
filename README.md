# Les astéroïdes

## Présentation

Sur le site du Jet Propulsion Laboratory, il est possible de récupérer des informations sur de nombreux asteroïdes, via l'application accessible en ligne ci-dessous.

https://ssd.jpl.nasa.gov/tools/sbdb_lookup.html#/

Les trajectoires de ces asteroïdes en provenance de ce site peuvent par ailleurs être téléchargées directement grâce au programme "asteroide.py" mis en ligne.

Le "notebook" : [ici](asteroide.ipynb)

Le programme "asteroide.py" : [ici](asteroide.py)

C'est l'occasion de revenir sur quelques astéroides qui ont fait l'actualité, notamment "2024 PT5", qualifié de "deuxième Lune temporairede la Terre", du fait de sa trajectoire assez particulière. Il s'agit en fait d'un asteroïde de la taille d'un gros bus, un beau caillou tout de même mais qui ne méritait pas le qualificatif de "Lune".

Plus anciennement, il y avait "Apophis" qui nous frôlera le 13 avril 2029.

Vidéo "Apophis" :

https://www.youtube.com/watch?v=l3eGlIXXz8Y

Et enfn, "Cruithne" qui suit une orbite "presque synchrone" avec la Terre, dite "en fer à cheval". C'est une situation particulière, associée aux points de Lagrange du système Soleil-Terre.

Vidéo "Les points de Lagrange": 

https://www.youtube.com/watch?v=KUQBsQVXbxE

On retrouvera aussi ces trajectoires avec des simulations numériques fondées sur la mécanique céleste.

## Description

### Introduction

Avec les programmes proposés, on pourra visualiser la trajectoire de différents astéroïdes, après avoir récupéré les éphémérides sur le site du JPL. Une simulation numérique de la trajectoire est ensuite envisagée dans le référentiel ICRS, avec la présence du Soleil et des huit planètes (et la Lune).

Les étapes décrites sont suivies avec les astéroïdes mentionnés ci-dessous.

- 2024 PT5
- Apophis
- Cruithne

### Etapes suivies

On génère le fichier avec la liste des astres qui nous concerneront.
```python
planetes(Planetes='Planetes.csv',Ligne = '2024PT5,0,km^3 s^-2,2024 PT5')
```

Ce fichier a le contenu visualisé ci-dessous.
```python
        Nom            mu      units        id
0    Soleil  1.327124e+11  km^3 s^-2        10
1   Mercure  2.203187e+04  km^3 s^-2         1
2     Venus  3.248586e+05  km^3 s^-2         2
3     Terre  3.986004e+05  km^3 s^-2       399
4      Lune  4.902800e+03  km^3 s^-2       301
5      Mars  4.282838e+04  km^3 s^-2         4
6   Jupiter  1.267128e+08  km^3 s^-2         5
7   Saturne  3.794058e+07  km^3 s^-2         6
8    Uranus  5.794556e+06  km^3 s^-2         7
9   Neptune  6.836527e+06  km^3 s^-2         8
10  2024PT5  0.000000e+00  km^3 s^-2  2024 PT5
```

Les éphémérides des astres présents dans la liste sont sauvegardées dans des fichiers CSV.
```python
convert_req_jpl_to_csv(Planetes='Planetes.csv',debut='2024-06-01 00:00:00',fin='2025-06-01 00:00:00',pas='24 h')
```

On peut alors visualiser la distance à laquelle l'astéroïde se trouve par rapport à un astre de référence, comme la Terre.
```python
trace2D(Astre1='Terre', Astre2='2024PT5')
```

La trajectoire de cet astéroïde peut aussi être visualisée.
```python
trace3D(Astre1='Terre', Astre2='2024PT5')
```

La simulation numérique est lancée.
```python
simu_aster_save(Planetes='Planetes.csv',methode='rk8')
```

La comparaison entre la trajectoire donnée par les éphémérides et la trajectoire simulée est alors possible.
```python
compare(Astre1_a='Terre',Astre2_a='2024PT5',Astre1_b='Terre',Astre2_b='2024PT5_simu')
```

---




