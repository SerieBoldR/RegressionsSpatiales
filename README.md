# Méthodes de régressions spatiales : un grand bol d’R

**Apparicio Philippe, Jérémy Gelb, Jean Dubé et Joan Carles Martori (2025). *Méthodes de régression spatiale : un grand bol d'R*. Université Laval et Université de Sherbrooke. fabriqueREL. Licence CC BY-SA.**

**Résumé :** Ce livre vise à décrire une panoplie de méthodes de régression spatiale avec le logiciel ouvert R. La philosophie de ce livre est de donner toutes les clés de compréhension et de mise en œuvre des méthodes abordées dans R. La présentation des méthodes est basée sur une approche compréhensive et intuitive plutôt que mathématique, sans pour autant négliger la rigueur statistique.

Le livre a d'ailleurs été écrit intégralement dans R avec [Quarto](https://quarto.org/).

## Structure du livre

Le manuel est structuré en six parties.

**Partie 1. Notions de base**

Dans cette première partie, nous présentons les jeux de données utilisés pour mettre en œuvre les différentes méthodes de régression spatiale présentées dans le livre. Nous discutons aussi de plusieurs notions fondamentales et méthodes qu'il importe de bien maîtriser avant d'aborder les chapitres suivants consacrés aux méthodes de régression spatiale, notamment l'autocorrélation spatiale, la notion de variable spatialement décalée et la régression linéaire multiple. Nous conclurons ce chapitre en exposant les raisons qui justifient l’utilisation de différentes formes de régressions spatiales pour modéliser des données spatiales ou spatiotemporelles.

**Partie 2. Spécification de la structure de covariance spatiale**

Dans cette seconde partie, les types de régression spatiale retenus introduisent l’espace en spécifiant directement la structure de covariance induite par l’autocorrélation spatiale dans les matrices de covariance de distribution normale. La méthode des moindres carrés généralisés (GLS) permet d’étendre le modèle des moindres carrés ordinaires (MCO) pour tenir compte de cette structure de covariance entre les observations. Cette approche est similaire à l’introduction d’effets aléatoires, ce qui a notamment conduit à la construction de GLMM (modèles linéaires généralisés à effets mixtes) introduisant des effets aléatoires distribués normalement avec une matrice de covariance structurée spatialement. Finalement, nous décrirons le modèle d'autorégression conditionnelle (*conditional autoregressive model*, CAR) qui s'applique à une variable dépendante continue dont la valeur pour une entité spatiale dépend de celles des entités spatiales proches ou voisines.

**Partie 3. Modèles d'économétrie spatiale**

Cette troisième partie comprend trois chapitres qui sont consacrés aux modèles d'économétrie spatiale qui vise à modéliser la **dépendance spatiale**. D'emblée, nous décrivons les principaux modèles spatiaux autorégressifs pour une variable dépendante continue qui permettent d'introduire l'autocorrélation spatiale sur les variables indépendantes (modèle SLX), la variable dépendante (SAR), le terme d'erreur (SEM), à la fois la variable dépendante et les variables indépendantes (SDM) et à la fois les variables indépendantes et le terme d’erreur (SDEM). Puis, nous abordons les modèles probit spatiaux pour modéliser une variable qualitative dichotomique (binaire). Finalement, nous décrivons d'autres extensions des modèles autorégressifs, soit les modèles spatiaux en panel qui permettent de modéliser des données spatiales longitudinales.

**Partie 4. Variable latente spatiale : lissage et filtrage spatial**

Dans cette troisième partie, les modèles retenus ont la particularité d’ajouter un terme spatialement structuré dans leur équation de régression. Ce terme spatial est construit à partir d’un ensemble de fonctions de base multipliées par des coefficients. Ces fonctions de base suivent des patrons géographiques, ce qui leur permet de capturer une variable latente spatialement autocorrélée qui autrement aurait fini dans les résidus. Le premier chapitre de cette partie est consacré aux modèles généralisés additifs (GAM) qui permettent d'introduire l’espace de deux manières différentes : avec une *spline* bivariée construite à partir des coordonnées géographiques (x, y) pour capturer les variations continues dans l'espace; avec un lissage par champ aléatoire de Markov (*Markov Random Field* – MRF) pour modéliser la dépendance spatiale entre les unités spatiales voisines. Dans le second, nous abordons les modèles linéaires généralisés avec des vecteurs spatiaux (*Spatial Eigenvector Generalized Linear Models*, SEVM). Ces modèles SEVM ajoutent des variables construites à partir de la décomposition de la matrice de pondération spatiale en vecteurs propres (*Moran Eigenvectors*). Ces deux types de modèles (GAM et SEVM) se ressemblent énormément dans leur conceptualisation de l’espace.

**Partie 5. Régressions spatiales et hétérogénéité spatiale**

Dans cette cinquième partie, nous abordons plusieurs types de régression spatiale qui permettent de faire varier les coefficients de régression dans l'espace. Premièrement, les régressions géographiquement pondérées (*Geographically weighted regression* - GWR) permettent d'explorer et de visualiser l'**hétérogénéité spatiale**, soit l'instabilité des relations entre la variable dépendante et les variables indépendantes. Premièrement, nous décrirons les formes dites classiques de la GWR qui s'appliquent à des variables dépendantes continues, dichotomiques (logistique) et de comptage (Poisson). Puis, nous abordons des extensions de la GWR, particulièrement la GWR mixte, la GWR multiéchelle et les GWR mixtes avec des variables spatialement décalées (MGWR-SAR). Finalement, nous verrons comment il est possible d'introduire des coefficients variant spatialement avec des modèles GAM (modèles généralisés additifs) et GLMM (modèles linéaires généralisés à effets mixtes).

**Partie 6. Conclusions**

Cette dernière partie regroupe les exercices corrigés, une conclusion générale et la bibliographie.
