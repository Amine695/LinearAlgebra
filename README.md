# Linear Algebra project

C project with the implementations of several algebraic algorithmic methods such as : <br/>

- The LU decomposition for a square matrix <br/>
- Solving a linear system using the LU decomposition <br/>
- The "naive" matrix inversion using the LU decomposition <br/>
- The inverse matrix using the Strassen's inversion algorithm with naive product <br/>
- Compute the inverse matrix but with the Strassen's multiplication algorithm <br/>
- A benchmark to compute time measurements over the different algorithms. 


For more details, please see the doxygen documentation opening the **index.html** file in a browser.


## Compilation
Un makefile est disponible pour la compilation du projet. Pour l'utiliser, executer la commande `make` dans le dossier principal. Un dossier **build** sera crée qui contiendra les fichiers objets, ainsi qu'un dossier **benchmarks** pour l'exportation des données.

## Execution 
### Algorithmes
Pour exécuter les algorithmes demandés, `./project --n --p` avec **n** la taille de la matrice (préferez une puissance de 2) et **p** un nombre premier.

### Benchmark
Pour exécuter les benchmarks, `./project benchmarks n` avec **n** la taille limite en puissance de 2 donc par exemple, n = 11 fait tourner les algorithmes sur des tailles allant de $2^1$ à $2^{11}$.
Une fois exécuter, on vous proposera de choisir la comparaison que vous souhaitez effectuer.<br/>
A la fin du benchmark, un **plot** de vos résultat sera généré automatiquement.<br/> Toutefois, si rien n'apparait ou qu'une erreur s'affiche, utiliser la commande `make plot` directement.<br/>
Compiler sous **Linux**. Si vous utilisez **python** au lieu de **python3** dans le shell, veuillez changer le chemin de l'interpréteur python en haut de `plot.py`. 

Authors: <br/>
Amine Berbagui <br/>
Ghassen Hachani