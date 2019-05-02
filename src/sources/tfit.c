#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

// Méthode de Transfert-Fusion itérée
// Lexique:
//   N = nbVertices
//   K = nbClasses
//   A = adjacency matrix
//   B = modularity matrix
//   P = partition

int Transfert(double** B, int N, int* P, int* K)
{
  // Calcul de la contribution de chaque élément à chaque classe
  double** Var = malloc(N*sizeof(double*));
  int* Clas = malloc(N*sizeof(int));
  for (int i=0; i<N; i++)
  {
    Var[i] = calloc(*K, sizeof(double));
    for (int j=0; j<N; j++)
      Var[i][P[j]] += B[i][j];
    // meilleure affectation dans Clas[i]
    double VarMax = -1.;
    int kk = 0;
    for (int k=0; k < *K; k++)
    {
      if (Var[i][k] > VarMax)
      {
        VarMax = Var[i][k];
        kk = k;
      }
    }
    if (VarMax > 0.)
      Clas[i] = kk;
    else
      Clas[i] = -1;
  }

  // fflag=1 si au moins 1 transfert dans cet appel de la procédure
  int NbTrans=0, OldC, NewC;
  while (1)
  {
    double gainmax = 0.;
    int ii = -1;
    // cherche le meilleur transfert
    for (int i=0; i<N; i++)
    {
      double gain = (Clas[i] >= 0 ? Var[i][Clas[i]]-Var[i][P[i]] : 0.);
      if (gain > gainmax)
      {
        gainmax = gain;
        ii = i;
      }
    }
    if (ii == -1) // plus rien à gagner
      break;
    NbTrans++; // on deplace ii ; mise a jour de Var
    int OldC = P[ii];
    int NewC = Clas[ii];

    for (int j=0; j<N; j++)
      Var[j][OldC] -= B[j][ii];
    if (NewC < 0)
    {
      (*K)++;
      NewC = *K-1;
    }
    P[ii] = NewC;
    for (int j=0; j<N; j++)
    {
      Var[j][NewC] += B[j][ii];
      if (Clas[j] == OldC)
      {
        // on recalcule la meilleure affectation
        double VarMax = -1.;
        int kk = 0;
        for (int k=1; k < *K; k++)
        {
          if (Var[j][k] > VarMax)
          {
            VarMax = Var[j][k];
            kk = k;
          }
        }
        if (VarMax>0.)
          Clas[j] = kk;
        else
          Clas[j] = -1;
      }
      else if (Var[j][NewC] > Var[j][Clas[j]])
        Clas[j] = NewC;
    }
  }

  for (int i=0; i<N; i++)
    free(Var[i]);
  free(Var);
  free(Clas);

  return NbTrans;
}

double Score(double** B, int N, int* P)
{
  double  Sc = 0.;
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      if (P[i] == P[j])
        Sc += B[i][j];
    }
  }
  return Sc;
}

void Around(double** B, int N, int* P, int* K)
{
  int* Q = malloc(N*sizeof(int)); // subdivision d'une classe
  for (int i=0; i<N; i++)
    Q[i] = P[i];
  int NbC = *K; //mémorise la meilleure partition
  double ScIni = Score(B, N, P);
  int NbEs = N;
  if (NbEs > 500) //TODO: explain this number
    NbEs = 500;
  int es = 0;
  while (es < NbEs)
  {
    for (int i=0; i<N; i++)
      P[i] = Q[i]; //on repart de la meilleure
    *K = NbC;
    if (*K == 1)
      return;
    int k = 0;
    // Nb de swapping diminue avec les essais:
    int kmax = 1 + (1.*rand()/RAND_MAX)*N / *K;
    while (k < kmax)
    {
      int i = (1.0*rand()/RAND_MAX)*N; //i indice au hasard entre 0 et N-1
      int j = (1.0*rand()/RAND_MAX)*N; //j indice au hasard entre 0 et N-1
      if (P[i] != P[j])
      {
        int ij = P[i];
        P[i] = P[j];
        P[j] = ij;
        k++;
      }
    }

    Transfert(B, N, P, K); //NOTE: unused result (NbTrans)
    double Sc = Score(B, N, P);
    if (Sc > ScIni)
    {
      for (int i=0; i<N; i++)
        Q[i] = P[i];
      ScIni = Sc;
      NbC = *K;
      es = 0;
    }
    else
      es++;
  }
  *K = NbC;
  for (int i=0; i<N; i++)
    P[i] = Q[i];

  free(Q);
  return;
}

int Louv1(double** B, int N, int* P, int* K)
{
  // Calcul de la contribution de chaque element a chaque classe
  double** Var = malloc(N*sizeof(double*));
  for (int i=0; i<N; i++)
  {
    Var[i] = calloc(N, sizeof(double)); //N, because K could increase
    for (int j=0; j<N; j++)
      Var[i][P[j]] += B[i][j];
  }

  // fflag = 1 si au moins 1 transfert dans cet appel de la procédure
  int flag = 1;
  int fflag = 0;
  while (flag > 0)
  {
    flag = 0; //transfert dans cette boucle
    for (int i=0; i<N; i++)
    {
      int OldC = P[i];
      double VarMax = Var[i][OldC];
      int NewC = -1;
      for (int k=0; k < *K; k++)
      {
        if (Var[i][k] > VarMax)
        {
          VarMax = Var[i][k];
          NewC = k;
        }
      }
      if (NewC < 0 && VarMax >= 0.)
        continue;
      flag = 1;
      fflag = 1;
      for (int j=0; j<N; j++)
        Var[j][OldC] -= B[j][i];
      if (VarMax < 0.)
      {
        (*K)++;
        NewC = *K - 1;
        for (int j=0; j<N; j++)
          Var[j][NewC] = 0.;
      }
      P[i] = NewC;
      for (int j=0; j<N; j++)
        Var[j][NewC] += B[j][i];
    }
  }

  for (int i=0; i<N; i++)
    free(Var[i]);
  free(Var);

  return fflag;
}

// Renumérote les classes
void Renum(double** B, int N, int* P, int* K)
{
  // Nb. d'elements dans les classes courantes:
  int* Kard = calloc(*K, sizeof(int));
  for (int i=0; i<N; i++)
    Kard[P[i]]++;
  int kk = 0;
  for (int k=0; k < *K; k++)
  {
    if (Kard[k] == 0)
      continue;
    kk++;
    for (int i=0; i<N; i++)
    {
      if (P[i] == k)
        P[i] = kk; //kk <= k, so it's OK
    }
  }
  *K = kk;
}

int Louv2(double** B, int N, int* P, int* K)
{
  // Poids des connections entre classes:
  double** W = malloc(*K * sizeof(double*));
  for (int k=0; k < *K; k++)
    W[k] = calloc(*K, sizeof(double));
  for (int i=0; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      W[P[i]][P[j]] += B[i][j];
      W[P[j]][P[i]] += B[i][j];
    }
  }

  double** Var = malloc(*K * sizeof(double*));
  int* Clas = malloc(*K * sizeof(int));
  for (int k=0; k < *K; k++)
  {
    Var[k] = malloc(N*sizeof(double)); //N, because K could increase
    Clas[k] = k;
    for (int kk=0; kk < *K; kk++)
      Var[k][kk] = W[k][kk];
  }

  // on fusionne des qu'il y a une connection > 0 entre classes
  int flag = 1;
  int fflag = 0;
  while (flag > 0)
  {
    flag = 0;
    for (int i=0; i < *K; i++)
    {
      int OldC = Clas[i];
      double VarMax = Var[i][OldC];
      int NewC = -1;
      for (int k=0; k < *K; k++)
      {
        if (Var[i][k] > VarMax)
        {
          VarMax = Var[i][k];
          NewC = k;
        }
      }
      if (VarMax < 0.)
      {
        (*K)++;
        NewC = *K - 1;
        for (int k=1; k < *K; k++)
        {
          Var[k][NewC] = 0.;
          Var[NewC][k] = W[k][i];
        }
        Var[NewC][OldC] = Var[i][OldC];
      }
      if (NewC < 0)
        continue;
      flag = 1;
      fflag = 1;
      Clas[i] = NewC;
      // on deplace i de OldC a NewC ; mise a jour de Var
      for (int j=0; j < *K; j++)
        Var[j][OldC] -= W[j][i];
      for (int j=0; j < *K; j++)
        Var[j][NewC] += W[j][i];
    }
  }

  if (fflag)
  {
    for (int i=0; i<N; i++)
      P[i] = Clas[P[i]];
  }

  for (int i=0; i < *K; i++)
  {
    free(Var[i]);
    free(W[i]);
  }
  free(Var);
  free(W);
  free(Clas);

  return fflag;
}

int TFit(double** B, int N, int* P)
{
  // Initialize partition:
  for (int i=0; i<N; i++)
    P[i] = i;
  int K = N;

  int NbPas = 0; //NOTE: for debug
  while (1)
  {
    NbPas++;
    if (!Louv1(B, N, P, &K))
      break;
    Renum(B, N, P, &K);
    // Y a-t-il des connections > 0 entre classes ?
    if (!Louv2(B, N, P, &K))
      break;
    Renum(B, N, P, &K);
  }

  //return NbPas;
  return K; //TODO
}

// Calcule la matrice des pondérations des paires
double** MatrixMod(double alpha, int N, int** A)
{
  double** B = malloc(N * sizeof(double*)); //matrice des modularités
  for (int i=0; i<N; i++)
    B[i] = malloc (N * sizeof(double));

  // On évalue les sommes en chaque sommet
  double SumMax = 0.;
  // Somme des poids des aretes en chaque sommet:
  double* Sum = malloc(N * sizeof(double));
  for (int i=0; i<N; i++)
  {
    Sum[i] = 0.;
    B[i][i] = 0.;
    for (int j=0; j<N; j++)
      Sum[i] += A[i][j];
    SumMax += Sum[i];
  }

  //  Matrice B
  double ModTot=0.;
  double ModMax=0.;
  for (int i=1; i<N; i++)
  {
    for (int j=0; j<i; j++)
    {
      B[i][j] = alpha*SumMax*A[i][j] - Sum[i]*Sum[j]/alpha; //ma formule
      // B[i][j]=(1+alpha)*SumMax*A[i][j]-Sum[i]*Sum[j]-alpha*SumMax; //celle de Fred
      // B[i][j]=(1+alpha)*SumMax*A[i][j]-alpha*SumMax; //new Fred
      B[j][i] = B[i][j];
      ModTot += B[i][j];
      if (A[i][j])
        ModMax += B[i][j];
      // La partie supérieure droite marquera à 1 les paires d'éléments
      // réunis dans au moins une classe.
    }
  }
  free(Sum);
  return B;
}

// Calcule une partition qui optimise un critère de modularité
// Algorithme de Transfert-Fusion itéré
void tfit_core(int** A, int N, int* P)
{
  double** B = MatrixMod(1.0, N, A);
  int K = TFit(B, N, P); //NOTE: unused result NbPas (for debug)
  Around(B, N, P, &K);
  //ClasOut(g, p);

  // Release memory
  for (int i=0; i<N; i++)
    free(B[i]);
  free(B);
}
