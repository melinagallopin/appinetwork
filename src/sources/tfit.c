#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

// Méthode de Transfert-Fusion itérée

typedef struct Graph
{
  char** Et; //étiquettes
  char** T; //graph table (adjacence)
  short* Dg; //degrés des sommets
  int Na; //nombre d'arêtes
  int N; //nombre de sommets
} Graph;

typedef struct Partition
{
  short *Clas; //les classes des sommets
  short *Part; //la partition
  int NbClas; //nombre de classes
} Partition;

Partition* AllocPartition(int N)
{
  Partition* p;
  short* Clas = malloc((N+1) * sizeof(short)); //classe de chaque element
  short* Part = malloc(N * sizeof(short)); //partition courante
  p->Clas = Clas;
  p->Part = Part;
  p->NbClas = N;
  return p;
}

// Lit le graphe, dans la table binaire T ; Calcule les degrés Dg[]
// TODO: à faire en R
Graph* LecGraph(const char* X)
{
  char Ch1[64],Ch2[64],OldCh[64]="";

  char FichE[1024]="Data/";
  strcpy(FichE, X);
  printf("%s\n",FichE);

  char aqw[512];
  getcwd(aqw, 512);
  printf("%s\n",aqw);

  FILE* FichCar = fopen(FichE,"r");
  assert(FichCar != NULL);

  // Première passe : On compte les sommets et on les range (ordre quelconque)
  int Na;
  fscanf(FichCar,"%d",&Na);
  char** Et_tmp = malloc(Na * sizeof(char*)); //Na >> N, but we don't know N
  strcpy(OldCh,"");
  int flag;
  int N = 0;
  for (int i=0; i<Na; i++)
  {
    fscanf(FichCar,"%s",Ch1);
    fscanf(FichCar,"%s",Ch2);
    if (strcmp(Ch1,OldCh) != 0) // si ce n'est pas le meme que precedemment
    {
      flag=1;
      for (int j=0; j<N; j++)
      {
        if (strcmp(Ch1,Et_tmp[j]) == 0)
        {
          flag=0;
          strcpy(OldCh,Ch1);
          break;
        }
      }
      if (flag)
      {
        Et_tmp[N] = malloc(strlen(Ch1)+1);
        strcpy(Et_tmp[N],Ch1);
        N++;
      }
    }
    flag=1;
    for (int j=0; j<N; j++)
    {
      if (strcmp(Ch2,Et_tmp[j]) == 0)
      {
        flag=0;
        break;
      }
    }
    if (flag)
    {
      Et_tmp[N] = malloc(strlen(Ch2)+1);
      strcpy(Et_tmp[N],Ch2);
      N++;
    }
  }
  fclose(FichCar);
  char** Et = malloc(N*sizeof(char*));
  for (int i=0; i<N; i++)
  {
    Et[i] = malloc(strlen(Et_tmp[i])+1);
    strcpy(Et[i], Et_tmp[i]);
    free(Et_tmp[i]);
  }
  free(Et_tmp);
  
  // seconde passe : on établit la table T d'adjacence
  char** T = malloc((N) * sizeof(char *));
  short* Dg = malloc((N) * sizeof(short));  // Degré de chaque sommet
  for (int i=0; i<N; i++)
  {
    T[i] = calloc ( N, sizeof(int) );
    Dg[i]=0;
  }
  
  FichCar = fopen(FichE,"r");
  fscanf(FichCar,"%d",&Na);
  strcpy(OldCh,"");
  int ii=-1;
  for (int i=0; i<Na; i++)
  {
    fscanf(FichCar,"%s",Ch1);
    fscanf(FichCar,"%s",Ch2);
    if (strcmp(Ch1,Ch2) == 0) continue;
    if (strcmp(Ch1,OldCh) != 0) // si ce n'est pas le meme que precedemment
    {
      for (int j=0; j<N; j++)
      {
        if (strcmp(Ch1,Et[j]) == 0)
        {
          ii=j;
          strcpy(OldCh,Ch1);
          break;
        }
      }
    }
    for (int j=0; j<N; j++)
    {
      if (strcmp(Ch2,Et[j]) == 0)
      {
        T[ii][j]=1;
        T[j][ii]=1; /* Le graphe est symétrisé */
        break;
      }
    }
  }
  fclose(FichCar);
  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      if (T[i][j])
      {
        Dg[i]++;
        Dg[j]++;
      }
    }
  }
  Graph* g = malloc(sizeof(Graph));
  g->Et = Et;
  g->T = T;
  g->Dg = Dg;
  g->Na = Na;
  g->N = N;
  return g;
}

// Calcule la matrice des pondérations des paires
float*** MatrixMod(float alpha, const Graph* g)
{
  int N = g->N;
  float** A = malloc(N * sizeof(float *)); //graphe pondéré
  float** B = malloc(N * sizeof(float *)); //matrice des modularités
  float** Var = malloc(N * sizeof(float *)); //variation en cas de fusion
  for (int i=0; i<N; i++)
  {
    A[i] = malloc (N * sizeof(float));
    B[i] = malloc (N * sizeof(float));
    Var[i] = malloc ((N+1) * sizeof(float));
    Var[i][i] = 0.;
    A[i][i] = 0.;
  }
  // Le graphe donné ponderé par 1 (TODO: ???)
  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      float val = (g->T[i][j] ? 1. : 0.);
      A[i][j]=val;
      A[j][i]=val;
    }
  }

  // On évalue les sommes en chaque sommet
  double SumMax=0.;
  // Somme des poids des aretes en chaque sommet:
  float* Sum = malloc(N * sizeof(float));
  for (int i=0; i<N; i++)
  {
    Sum[i]=0.;
    B[i][i]=0.;
    for (int j=0; j<N; j++)
    {
      if (A[i][j]>0.)
        Sum[i]=Sum[i]+A[i][j];
    }
    SumMax += Sum[i];
  }

  //  Matrice B
  double ModTot=0.;
  double ModMax=0.;
  SumMax /= 1.; //TODO: ???
  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      B[i][j]=alpha*SumMax*A[i][j]-Sum[i]*Sum[j]/alpha; //ma formule
      // B[i][j]=(1+alpha)*SumMax*A[i][j]-Sum[i]*Sum[j]-alpha*SumMax; //celle de Fred
      // B[i][j]=(1+alpha)*SumMax*A[i][j]-alpha*SumMax; //new Fred
      B[j][i]=B[i][j];
      ModTot += B[i][j];
      if (A[i][j] > 0.)
        ModMax += B[i][j];
      // La partie supérieure droite marquera à 1 les paires d'éléments
      // réunis dans au moins une classe.
      A[j][i] = 0.;
    }
  }
  free(Sum);
  float*** res = malloc(3*sizeof(float**));
  res[0] = A;
  res[1] = B;
  res[2] = Var;
  return res;
}

int Louv1(const float** B, float** Var, int N, Partition* p)
{
  // Calcul de la contribution de chaque element a chaque classe
  for (int i=0; i<N; i++)
  {
    for (int k=1; k<=p->NbClas; k++)
      Var[i][k] = 0.;
    for (int j=0; j<N; j++)
    {
      int k = p->Part[j];
      Var[i][k] += B[i][j];
    }
  }

  // fflag = 1 si au moins 1 transfert dans cet appel de la procédure
  int flag = 1;
  int fflag = 0;
  while (flag>0)
  {
    flag = 0;  //transfert dans cette boucle
    for (int i=0; i<N; i++)
    {
      int OldC = p->Part[i];
      float VarMax = Var[i][OldC];
      int NewC = 0;
      for (int k=1; k <= p->NbClas; k++)
      {
        if (Var[i][k] > VarMax)
        {
          VarMax = Var[i][k];
          NewC = k;
        }
      }
      if (NewC <= 0 && VarMax >= 0.)
        continue;
      flag = 1;
      fflag = 1;
      for (int j=0; j<N; j++)
        Var[j][OldC] -= B[j][i];
      if (VarMax < 0.)
      {
        p->NbClas++;
        NewC = p->NbClas;
        for (int j=0; j<N; j++)
          Var[j][NewC] = 0.;
      }
      p->Part[i] = NewC;
      for (int j=0; j<N; j++)
        Var[j][NewC] += B[j][i];
    }
  }
  return fflag;
}

// Renumérote les classes
void Renum(const float** B, int N, Partition* p)
{
  // Nb. d'elements dans les classes courantes:
  short* Kard = calloc((p->NbClas+1), sizeof(short));
  for (int i=0; i<N; i++)
    Kard[p->Part[i]]++;
  int kk = 0;
  for (int k=1; k <= p->NbClas; k++)
  {
    if (Kard[k]==0)
      continue;
    else
      kk++;
    int card=0;
    for (int i=0; i<N; i++)
    {
      if (p->Part[i] == k)
      {
        p->Clas[card] = i;
        card++;
      }
    }
    for (int i=0; i<card; i++)
      p->Part[p->Clas[i]] = kk;
  }
  p->NbClas=kk;

  // TODO: ancien contenu "utile" de ClasOut()
  // Renumérotation des classes de 1 à NbClas
  int flag = 1;
  while (flag)
  {
    flag=0;
    for (k=1; k<=p->NbClas; k++)
    {
      if (Kard[k]>0)
        continue;
      for (i=0; i<g->N; i++)
      {
        if (p->Part[i]==p->NbClas)
          p->Part[i]=k;
      }
      Kard[k]=Kard[p->NbClas];
      p->NbClas--;
      flag=1;
    }
  }

  TotAint=0; // Modularite du graphe non pondere;
  sing=0; //nb. aretes intraclasse;
  for (k=1; k<=p->NbClas; k++)
  {
    printf("Class %3d Nb. of elements %d\n",k,Kard[k]);
    card=0;
    DgInt=0; //somme des degrés des sommets de la classe
    for (i=0; i<g->N; i++)
    {
      if (p->Part[i]==k)
      {
        p->Clas[card]=i;
        card++;
        DgInt += g->Dg[i];
        printf("%s ",g->Et[i]);
      }
    }
}

int Louv2(float** A, const float** B, float** Var, int N, Partition* p)
{
  // poids des connections entre classes
  for (int k=0; k <= p->NbClas; k++)
  {
    for (int kk=0; kk <= p->NbClas; kk++)
      A[k][kk] = 0.;
  }
  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      A[p->Part[i]][p->Part[j]] += B[i][j];
      A[p->Part[j]][p->Part[i]] += B[i][j];
    }
  }

  for (int k=1; k<=p->NbClas; k++)
  {
    p->Clas[k] = k;
    A[k][k] = 0.;
    for (int kk=1; kk <= p->NbClas; kk++)
      Var[k][kk] = A[k][kk];
  }

  // on fusionne des qu'il y a une connection > 0 entre classes
  int flag = 1;
  int fflag = 0;
  while(flag>0)
  {
    flag=0;
    for (int i=1; i<=p->NbClas; i++)
    {
      int OldC = p->Clas[i];
      float VarMax = Var[i][OldC];
      int NewC = 0;
      for (int k=1; k<=p->NbClas; k++)
      {
        if (Var[i][k] > VarMax)
        {
          VarMax = Var[i][k];
          NewC = k;
        }
      }
      if (VarMax < 0.)
      {
        p->NbClas++;
        NewC = p->NbClas;
        for (int k=1; k<p->NbClas; k++)
        {
          Var[k][p->NbClas] = 0.;
          Var[p->NbClas][k] = A[k][i];
        }
        Var[p->NbClas][OldC] = Var[i][OldC];
      }
      if (NewC <= 0)
        continue;
      flag = 1;
      fflag = 1;
      p->Clas[i] = NewC;
      // on deplace i de OldC a NewC ; mise a jour de Var
      for (int j=1; j <= p->NbClas; j++)
        Var[j][OldC] -= A[j][i];
      for (int j=1; j <= p->NbClas; j++)
        Var[j][NewC] += A[j][i];
    }
  }
  return fflag;
}

int TFit(float** A, const float** B, float** Var, int N, Partition* p)
{
  // Initialize partition:
  for (int i=0; i<N; i++)
    p->Part[i] = i+1;
  p->NbClas = N;

  int NbPas=0;
  while (1)
  {
    NbPas++;
    int fflag = Louv1(B, Var, N, p);
    if (!fflag)
      break;
    Renum(B, N, p);
    // Y a-t-il des connections > 0 entre classes ?
    fflag = Louv2(A, B, Var, N, p);
    if (!fflag)
      break;
    for (int i=0; i<N; i++)
      p->Part[i]=p->Clas[p->Part[i]];
    Renum(B, N, p);
  }
  return NbPas;
}

int Transfert(const float** B, float** Var, int N, Partition* p)
{
  // Calcul de la contribution de chaque élément à chaque classe
  for (int i=0; i<N; i++)
  {
    for (int k=1; k <= p->NbClas; k++)
      Var[i][k] = 0.;
    for (int j=0; j<N; j++)
    {
      int k = p->Part[j];
      Var[i][k] += B[i][j];
    }
    // meilleure affectation dans Clas[i]
    float VarMax = -1.;
    int kk = 0;
    for (int k=1; k <= p->NbClas; k++)
    {
      if (Var[i][k] > VarMax)
      {
        VarMax = Var[i][k];
        kk = k;
      }
    }
    if (VarMax > 0.)
      p->Clas[i] = kk;
    else
      p->Clas[i] = 0;
  }

  // fflag=1 si au moins 1 transfert dans cet appel de la procédure
  int flag = 1;
  int NbTrans=0, OldC, NewC;
  while (flag>0)
  {
    float gainmax=0.;
    int ii=-1;
    // cherche le meilleur transfert
    for (int i=0; i<N; i++)
    {
      float gain=Var[i][p->Clas[i]]-Var[i][p->Part[i]];
      if (gain>gainmax)
      {
        gainmax=gain;
        ii=i;
      }
    }
    if (ii==-1)
    {
      flag=0;
      continue; // plus rien à gagner
    }
    NbTrans++; // on deplace ii ; mise a jour de Var
    OldC=p->Part[ii];
    NewC=p->Clas[ii];
    
    for (int j=0; j<N; j++)
      Var[j][OldC] -= B[j][ii];
    if (NewC==0)
    {
      p->NbClas++;
      NewC=p->NbClas;
    }
    p->Part[ii]=NewC;
    for (int j=0; j<N; j++)
    {
      Var[j][NewC] += B[j][ii];
      if (p->Clas[j]==OldC) // on recalcule la meilleure affectation
      {
        float VarMax=-1.;
        int kk = 0;
        for (int k=1; k<=p->NbClas; k++)
        {
          if (Var[j][k]>VarMax)
          {
            VarMax=Var[j][k];
            kk=k;
          }
        }
        if (VarMax>0.)
          p->Clas[j]=kk;
        else
          p->Clas[j]=0;
      }
      else if (Var[j][NewC] > Var[j][p->Clas[j]])
        p->Clas[j]=NewC;
    }
  }
  return NbTrans;
}

double Score(const float** B, int N, short* Part)
{
  double  Sc=0.;
  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      if (Part[i]==Part[j])
        Sc += B[i][j];
    }
  }
  return Sc;
}

void Around(const float** B, float** Var, int N, Partition* p)
{
  short* Q = malloc(N * sizeof(short)); // subdivision d'une classe 
  for (int i=0; i<N; i++)
    Q[i] = p->Part[i];
  int NbC = p->NbClas; //mémorise la meilleure partition
  double ScIni = Score(B, N, p->Part);
  printf("Sc Ini %.0f : ",ScIni);
  int NbEs = N;
  if (NbEs > 500) //TODO: explain this number
    NbEs = 500;
  int es = 0;
  while (es < NbEs)
  {
    for (int i=0; i<N; i++)
      p->Part[i] = Q[i]; //on repart de la meilleure
    p->NbClas = NbC;
    if (p->NbClas == 1)
      return;
    int k=0;
    // Nb de swapping diminue avec les essais:
    int kmax = 1 + (1.*rand()/RAND_MAX)*N / p->NbClas;
    while (k < kmax)
    {
      int i = (1.0*rand()/RAND_MAX)*N; //i indice au hasard entre 0 et N-1
      int j = (1.0*rand()/RAND_MAX)*N; //j indice au hasard entre 0 et N-1
      if (p->Part[i] != p->Part[j])
      {
        int ij = p->Part[i];
        p->Part[i]=p->Part[j];
        p->Part[j]=ij;
        k++;
      }
    }

    Transfert(B, Var, N, p); //NOTE: unused result (NbTrans)
    double Sc = Score(B, N, p->Part);
    if (Sc > ScIni)
    {
      for (int i=0; i<N; i++)
        Q[i] = p->Part[i];
      ScIni = Sc;
      NbC = p->NbClas;
      es = 0;
    }
    else
      es++;
  }
  p->NbClas = NbC;
  for (int i=0; i<N; i++)
    p->Part[i] = Q[i];
  return;
}

void SaveClas(const char* FichS, const Graph* g, Partition* p)
{
  //printf("Minimum cardinality of selected classes (0=all) ");
  //scanf("%d",&Kmin); //TODO: could be a parameter
  int Kmin = 0;
  FILE* FichClas = fopen(FichS,"w");
  assert(FichClas != NULL);

  // Nb. d'elements dans les classes courantes:
  short* Kard = calloc(p->NbClasi+1, sizeof(short));
  for (int i=0; i < g->N; i++)
    Kard[p->Part[i]]++;
  for (int k=1; k <= p->NbClas; k++)
  {
    if (Kard[k] >= Kmin)
      kk++;
  }
  fprintf(FichClas," %d\n",kk);

  for (int k=1; k <= p->NbClas; k++)
  {
    if (Kard[k] < Kmin)
      continue;
    fprintf(FichClas," %d\n",Kard[k]);
    for (int i=0; i < g->N; i++)
    {
      if (p->Part[i] == k)
        fprintf(FichClas,"%s ", g->Et[i]);
    }
    fprintf(FichClas,"\n\n");
  }
  fprintf(FichClas,"\n\n");
  fclose(FichClas);
}

// calcule une partition qui optimise un critère de modularité
// Algorithme de Transfert-Fusion itéré
void tfit_core(int* adjacencyMatrix, int nbVertices, int* partition)
{
  // Read graph data (edges, labels...)
  const Graph* g = LecGraph(X);
  Partition* p = AllocPartition(g->N);

  float*** mmod = MatrixMod(1.0, g->N);
  float** A = mmod[0];
  const float** B = mmod[1];
  float** Var = mmod[2];
  TFit(A, B, Var, g->N, p); //NOTE: unused result NbPas
  Around(B, Var);
  ClasOut(g, p);
  SaveClas(out, g, p); //sauvegarde du fichier Clas

  // Release memory
  for (int i=0; i<g->N; i++)
  {
    free(g->Et[i]);
    free(g->T[i]);
  }
  free(g->Et);
  free(g->T);
  free(g->Dg);
  free(p->Clas);
  free(p->Part);
}
