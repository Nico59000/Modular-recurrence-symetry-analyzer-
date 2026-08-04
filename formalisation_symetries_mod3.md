# Symétries de demi-période pour les suites modulo 3

## 1. Objet

Cette note formalise les phénomènes décrits dans le dépôt
`59200/k-bonacci-pisano-finder` sous les expressions informelles
« périodes complémentaires à trois » et « auto-complémentaires lorsqu'elles
sont coupées en deux ».

Le cadre invariant est celui des mots périodiques sur le corps

\[
\mathbb F_3=\{0,1,2\},
\]

où l'involution additive est

\[
\nu(x)=-x\pmod 3,
\qquad
\nu(0)=0,\quad \nu(1)=2,\quad \nu(2)=1.
\]

L'expression « complémentaire à 3 » doit donc être remplacée par
**opposé additif modulo 3** ou **complément par négation modulo 3**.
La somme décimale des représentants n'est pas toujours égale à \(3\),
puisque \(0\) est envoyé sur \(0\).

---

## 2. Trois propriétés distinctes

Soit \(w=(w_0,\ldots,w_{T-1})\in\mathbb F_3^T\) un mot de période
primitive \(T\).

### 2.1 Paire de périodes opposées

Deux périodes \(w\) et \(v\) sont une paire opposée si

\[
v_n=-w_n\pmod 3
\]

pour tout \(n\). Elles peuvent appartenir à deux cycles distincts.

### 2.2 Antipériodicité de demi-période

Le mot \(w\) est **antipériodique à demi-période** si \(T=2h\) et

\[
w_{n+h}=-w_n\pmod 3
\]

pour tout \(n\). Il s'écrit alors

\[
w=H\Vert(-H),
\qquad
H=(w_0,\ldots,w_{h-1}).
\]

Exemples :

\[
01120221=0112\Vert0221,
\]

\[
011022=011\Vert022,
\]

\[
01220211=0122\Vert0211.
\]

### 2.3 Renversement de la seconde moitié

L'opération

\[
\mathcal R(a_0,\ldots,a_{h-1})
=(a_{h-1},\ldots,a_0)
\]

est indépendante de la négation. Dans l'exemple de Fibonacci,

\[
\mathcal R(0221)=1220.
\]

Ainsi, \(1220\) est le **mot antipodal renversé** de \(0112\), et non la
seconde moitié elle-même.

---

## 3. Théorème de l'unique décalage opposé

Soit \(w\) un mot primitif non nul de période \(T\). Supposons qu'il
existe un décalage \(r\) tel que

\[
w_{n+r}=-w_n
\]

pour tout \(n\). Alors \(T\) est pair et

\[
r\equiv \frac T2\pmod T.
\]

### Preuve

En appliquant deux fois le décalage,

\[
w_{n+2r}=w_n.
\]

La primitivité implique \(T\mid 2r\). Le décalage \(r\) n'est pas nul
modulo \(T\), car un mot non nul sur \(\mathbb F_3\) ne satisfait pas
\(w=-w\). L'élément non trivial d'ordre \(2\) du groupe cyclique
\(\mathbb Z/T\mathbb Z\) existe uniquement lorsque \(T\) est pair et vaut
\(T/2\).

### Conséquence

Une période primitive est soit :

1. envoyée par la négation sur un cycle distinct ;
2. soit invariante par négation, auquel cas elle est nécessairement
   antipériodique à demi-période.

Cette dichotomie corrige la confusion entre les couples

\[
00111201,\quad00222102
\]

et les périodes intrinsèquement antipériodiques comme \(011022\).

---

## 4. Conséquences combinatoires

Si

\[
w=H\Vert(-H),
\]

alors :

\[
\#\{n:w_n=1\}=\#\{n:w_n=2\},
\]

\[
\#\{n:w_n=0\}\equiv0\pmod2,
\]

et

\[
\sum_{n=0}^{T-1}w_n=0
\quad\text{dans }\mathbb F_3.
\]

Ces conditions sont nécessaires, mais non suffisantes.

En utilisant le relèvement centré

\[
\widetilde 0=0,\qquad
\widetilde 1=1,\qquad
\widetilde 2=-1,
\]

le polynôme générateur s'écrit

\[
W(X)
=
\sum_{n=0}^{2h-1}\widetilde w_nX^n
=
(1-X^h)
\sum_{n=0}^{h-1}\widetilde w_nX^n.
\]

---

## 5. Caractérisation de Fourier

Soit

\[
\widehat w(k)
=
\sum_{n=0}^{2h-1}
\widetilde w_n
e^{-2\pi i kn/(2h)}.
\]

Alors

\[
w_{n+h}=-w_n
\]

si et seulement si

\[
\widehat w(k)=0
\]

pour tout indice pair \(k\).

### Preuve directe

Lorsque \(w_{n+h}=-w_n\),

\[
\widehat w(k)
=
\left(1-(-1)^k\right)
\sum_{n=0}^{h-1}
\widetilde w_n e^{-2\pi i kn/(2h)}.
\]

Le facteur est nul pour \(k\) pair.

Réciproquement, si tous les modes pairs sont nuls, le mot appartient au
sous-espace engendré par les modes impairs. Le décalage de \(h\) agit sur
le mode \(k\) par \((-1)^k=-1\), donc il agit par \(-I\) sur le mot.

Cette caractérisation donne un diagnostic spectral exact, distinct d'une
simple comparaison visuelle des deux moitiés.

---

## 6. Suites linéaires modulo 3

Considérons une récurrence d'ordre \(k\)

\[
u_{n+k}
=
c_0u_n+c_1u_{n+1}+\cdots+c_{k-1}u_{n+k-1}
\pmod3.
\]

Son état est

\[
s_n=(u_n,\ldots,u_{n+k-1})^\mathsf T
\]

et son évolution est

\[
s_{n+1}=Ms_n,
\]

où \(M\) est la matrice compagnon.

### 6.1 Périodicité pure

Le déterminant de \(M\) vaut, au signe près, \(c_0\). Sur
\(\mathbb F_3\),

\[
M\text{ est inversible}
\iff
c_0\neq0.
\]

Dans ce cas, toute orbite d'état est purement périodique. Si \(c_0=0\),
une prépériode peut apparaître et doit être distinguée du cycle.

### 6.2 Négation des cycles

La linéarité donne

\[
M(-s)=-M(s).
\]

La négation envoie donc tout cycle sur un cycle de même longueur. Deux
cas seulement sont possibles :

- deux cycles distincts échangés par négation ;
- un même cycle invariant, nécessairement antipériodique à demi-période.

### 6.3 Critère sur une graine

Une graine \(s_0\) possède une antipériode \(h\) si

\[
M^hs_0=-s_0.
\]

Toute sortie linéaire de l'état satisfait alors

\[
u_{n+h}=-u_n.
\]

### 6.4 Critère uniforme

Toutes les graines possèdent la même antipériode \(h\) si et seulement si

\[
M^h=-I.
\]

Équivalemment, le polynôme minimal \(\mu_M\) vérifie

\[
\mu_M(X)\mid X^h+1
\quad\text{dans }\mathbb F_3[X].
\]

Il en résulte

\[
M^{2h}=I.
\]

---

## 7. Exemple de Fibonacci

Pour

\[
M=
\begin{pmatrix}
0&1\\
1&1
\end{pmatrix}
\quad\text{sur }\mathbb F_3,
\]

on calcule

\[
M^4=-I,
\qquad
M^8=I.
\]

La période issue de la graine \((0,1)\) est donc

\[
01120221
=
0112\Vert(-0112).
\]

Sa seconde moitié renversée est

\[
\mathcal R(0221)=1220.
\]

La coïncidence

\[
1220=L_{14}+F_{14}=2F_{15}
\]

est une identité d'encodage supplémentaire. Elle n'est pas une
conséquence générale de l'antipériodicité.

---

## 8. Généralisation aux suites arbitraires modulo 3

Aucune hypothèse de récurrence n'est nécessaire pour analyser un mot
périodique donné. Pour toute période primitive, on peut calculer :

- son cycle sous les rotations ;
- son opposé additif ;
- ses symétries affines
  \[
  w_{n+r}=aw_n+b;
  \]
- ses symétries de renversement
  \[
  w_{r-n}=aw_n+b;
  \]
- son éventuelle antipériodicité ;
- ses modes de Fourier pairs et impairs.

La théorie matricielle intervient uniquement lorsqu'une récurrence
linéaire déclarée est disponible.

---

## 9. Extension au module général

Sur \(\mathbb Z/m\mathbb Z\), l'involution naturelle reste

\[
x\longmapsto-x.
\]

La définition

\[
w_{n+h}=-w_n\pmod m
\]

reste valable. Pour \(m=2\), la négation est l'identité et la notion
d'antipériode se confond avec une période ordinaire. Pour \(m>2\), la
distinction demeure non triviale.

Le terme « complément à \(m\) » doit être évité comme définition
algébrique, car il dépend du choix des représentants entiers. Le terme
invariant est « opposé additif modulo \(m\) ».

---

## 10. Audit algorithmique

Une période d'une récurrence d'ordre \(k\) modulo \(m\) doit être
détectée sur l'espace fini des états

\[
(\mathbb Z/m\mathbb Z)^k,
\]

et non par recherche d'un bloc répété dans un préfixe de longueur
arbitraire.

Le script fourni utilise :

- le codage compact des états en base \(m\) ;
- l'algorithme de Brent pour une graine :
  \[
  O(\mu+\lambda)\text{ temps},\quad O(1)\text{ mémoire};
  \]
- une décomposition exacte du graphe fonctionnel pour toutes les graines :
  \[
  O(m^k)\text{ temps};
  \]
- la parallélisation entre familles de coefficients ;
- une classification séparée des cycles opposés et des cycles
  antipériodiques ;
- un contrôle matriciel de \(M^h=-I\).

Aucun noyau Mathematica ni wrapper externe n'est requis.

---

## 11. Résultats reproductibles du balayage

Le balayage des quinze familles nommées dans le dépôt retrouve notamment,
pour la famille tritetranacci, les classes suivantes :

- paires opposées distinctes de longueurs \(24\) et \(8\) ;
- périodes antipériodiques \(011022\), \(01220211\) et \(12\).

Un second balayage porte sur les \(119\) vecteurs de coefficients
binaires non nuls d'ordres \(2\) à \(6\) :

\[
\sum_{k=2}^{6}(2^k-1)=119.
\]

Les calculs produits par le script trouvent :

- \(13\) familles satisfaisant un critère global \(M^h=-I\) ;
- \(88\) familles possédant au moins un cycle antipériodique non trivial ;
- \(96\) familles possédant au moins une paire de cycles opposés distincts.

Ces nombres sont des résultats expérimentaux exhaustifs dans cette
fenêtre finie, reproductibles par les fichiers JSON et CSV fournis.
