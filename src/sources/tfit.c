#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <unistd.h>

#define MaxV 1500    // Nb. max de sommets
#define SupCar 21    // Nb. max de caracteres pour une etiquette

short  *Dg, *Clas, *Kard, *Part, *Q, *PartMax;
int      N=0, Na=0, NbClas, NbCar, CardMax=0;  // cardinal maximum d'une classe
char  Et[MaxV][SupCar], *string, **T;
float  *Sum, **A, **B, **Var, SumDg2=0., eps=.0001, alpha;
int    MCV=0, FuStyl=0;  // MVC=1 ssi recherche des affectations multiples ; FuStyl=0 si gain moyen, 1 si gain total

// Méthode de Transfert-Fusion itérée

static int compar(const void *e1, const void *e2)
{
  return strcmp((char*)e1,(char*)e2);
}

void EditGraph()
{
  int i,j;
  printf("Graphe\n");
  for (i=0; i<N; i++)
  {
    printf("%*s %3d : ",-NbCar,Et[i],Dg[i]);
    for (j=0; j<N; j++)
      printf("%1d ", T[i][j]);
    printf("\n");
  }
}

// Edite les classes a partir d'un vecteur de numeros de classe Part[]
// Calcule la modularite et le taux d'arêtes intraclasses
void ClasOut()
{
  int i, j1, j2, k, k1, k2, Aint, TotAint, DgInt, card, flag=1, sing;
  double ExpDen, Gain, SumGain=0., Dens, SumDens=0.;

  for (k=0; k<=NbClas; k++)
    Kard[k]=0;
  for (i=0; i<N; i++)
    Kard[Part[i]]++;

  // Renumérotation des classes de 1 à NbClas
  while (flag)
  {
    flag=0;
    for (k=1; k<=NbClas; k++)
    {
      if (Kard[k]>0)
        continue;
      for (i=0; i<N; i++)
      {
        if (Part[i]==NbClas)
          Part[i]=k;
      }
      Kard[k]=Kard[NbClas];
      NbClas--;
      flag=1;
    }
  }

  TotAint=0; // Modularite du graphe non pondere;
  sing=0; //nb. aretes intraclasse;
  for (k=1; k<=NbClas; k++)
  {
    printf("Class %3d Nb. of elements %d\n",k,Kard[k]);
    card=0;
    DgInt=0; //somme des degrés des sommets de la classe
    for (i=0; i<N; i++)
    {
      if (Part[i]==k)
      {
        Clas[card]=i;
        card++;
        DgInt += Dg[i];
        printf("%s ",Et[i]);
      }
    }
    ExpDen = 1.*DgInt / Na / 2;
    ExpDen = ExpDen*ExpDen;
    if (card<2)
    {
      printf("\n");
      sing++;
      continue;
    }
    Aint = 0;
    for (k1=1; k1<card; k1++)
    {
      j1=Clas[k1];
      for (k2=0; k2<k1; k2++)
      {
        j2=Clas[k2];
        if (T[j1][j2]>0) Aint++;
      }
    }
    Gain = 1.*Aint / Na - ExpDen;
    SumGain = SumGain + Gain;
    printf(" : mod %.3f\n",Gain);
    Dens = 2.*Aint/Kard[k]/(Kard[k]-1);
    SumDens = SumDens + Kard[k]*Dens;
    TotAint += Aint;
  }
  printf("\n");
  printf("Densite moy. %.3f  Modularity %.5f\n",SumDens/(N-sing),SumGain);
  printf("Nb. de singletons %d\n",sing);
  printf("Rate of intra class edges %.3f\n\n",1.*TotAint/Na);
}

// Lit le graphe, dans la table binaire T
// Calcule les degrés Dg[] et la somme des carrés des degrés SumDg2
void LecGraph(const char* X)
{
  int i, ii, j, jj, flag, NL;
  char MaxCar=0, Ch1[2*SupCar],Ch2[2*SupCar],OldCh[SupCar]="", car='N';
  float  seuil;

  char FichE[1024]="Data/";
  strcpy(FichE, X);

  // Non pondéré ou pondéré (TODO: this info should be in graph file)
  if (car=='N')
    strcat(FichE,".gr");
  else if (car=='O')
    strcat(FichE,".grp");
  printf("%s\n",FichE);
  char aqw[512];
  getcwd(aqw, 512);
  printf("%s\n",aqw);

  FILE* FichCar = fopen(FichE,"r");
  assert(FichCar != NULL);

  // Première passe : On compte les sommets et on les range par odre alphabetique
  fscanf(FichCar,"%d",&NL);
  strcpy(OldCh,""); NbCar=0;
  for (i=0; i<NL; i++)
  {
    fscanf(FichCar,"%s",Ch1);
    fscanf(FichCar,"%s",Ch2);
    if (car=='O')
      fscanf(FichCar,"%f",&seuil);
    if (strlen(Ch1)>NbCar)
      NbCar=strlen(Ch1);
    if (strlen(Ch2)>NbCar)
      NbCar=strlen(Ch2);
    if (NbCar>MaxCar)
      MaxCar=NbCar;
    if (MaxCar>=SupCar)
    {
      printf ("Labels are limited to %d characters\n",(SupCar-1));
      printf("Non correct edge : %5d   %s, %s\n",i,Ch1,Ch2);
      exit(0);
    }
    if (strcmp(Ch1,OldCh) != 0) // si ce n'est pas le meme que precedemment
    {
      flag=1;
      for (j=0; j<N; j++)
      {
        if (strcmp(Ch1,Et[j]) == 0)
        {
          flag=0;
          strcpy(OldCh,Ch1);
          break;
        }
      }
      if (flag)
      {
        strcpy(Et[N],Ch1);
        N++;
        if (N==MaxV)
        {
          printf("Nb. of vertices > %d  ",MaxV);
          getchar();
          return;
        }
      }
    }
    flag=1;
    for (j=0; j<N; j++)
    {
      if (strcmp(Ch2,Et[j]) == 0)
      {
        flag=0;
        break;
      }
    }
    if (flag)
    {
      strcpy(Et[N],Ch2);
      N++;
      if (N==MaxV)
      {
        printf("Nb. of vertices > %d  ",MaxV);
        getchar();
        exit(0);
      }
    }
  }
  fclose(FichCar);
  qsort(Et, N, SupCar, compar);
  
  // seconde passe : on établit la table T d'adjacence
  T = malloc((N) * sizeof(char *));
  Dg = malloc((N) * sizeof(short));  // Degré de chaque sommet
  assert (T != NULL && Dg != NULL);
  for (i=0; i<N; i++)
  {
    T[i] = malloc ( N * sizeof(char) );
    assert ( T[i] != NULL);
    Dg[i]=0;
    for (j=0; j<N; j++)
      T[i][j]=0;
  }
  
  FichCar = fopen(FichE,"r");
  fscanf(FichCar,"%d",&NL);
  strcpy(OldCh,""); ii=-1; Na=0;
  for (i=0; i<NL; i++)
  {
    fscanf(FichCar,"%s",Ch1);
    fscanf(FichCar,"%s",Ch2);
    if (car=='O')
      fscanf(FichCar,"%f",&seuil);
    if (strcmp(Ch1,Ch2) == 0) continue;
    if (strcmp(Ch1,OldCh) != 0) // si ce n'est pas le meme que precedemment
    {
      for (j=0; j<N; j++)
      {
        if (strcmp(Ch1,Et[j]) == 0)
        {
          ii=j;
          strcpy(OldCh,Ch1);
          break;
        }
      }
    }
    for (j=0; j<N; j++)
    {
      if (strcmp(Ch2,Et[j]) == 0)
      {
        jj=j;
        break;
      }
    }
    T[ii][jj]=1;
    T[jj][ii]=1; /* Le graphe est symetrisé */
  }
  fclose(FichCar);
  for (i=1; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      if (T[i][j]==1)
      {
        Dg[i]++;
        Dg[j]++;
        Na++;
      }
    }
  }
  for (i=0; i<N; i++)
    SumDg2=SumDg2+1.*Dg[i]*Dg[i];
  if (N<21)
    EditGraph();
}

// Calcule la matrice des pondérations des paires
void MatrixMod(float alpha)
{
  int i, j;
  double  SumMax, ModTot, ModMax;
  
  // On évalue les sommes en chaque sommet
  SumMax=0.;
  for (i=0; i<N; i++)
  {
    Sum[i]=0.;
    B[i][i]=0.;
    for (j=0; j<N; j++)
    {
      if (A[i][j]>0.)
        Sum[i]=Sum[i]+A[i][j];
    }
    SumMax=SumMax+Sum[i];
  }
  
  //  Matrice B
  ModTot=0.;
  ModMax=0.;
  SumMax=SumMax/1.;
  for (i=1; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      B[i][j]=alpha*SumMax*A[i][j]-Sum[i]*Sum[j]/alpha; // ma formule
      // B[i][j]=(1+alpha)*SumMax*A[i][j]-Sum[i]*Sum[j]-alpha*SumMax; // celle de Fred
      // B[i][j]=(1+alpha)*SumMax*A[i][j]-alpha*SumMax; // new Fred
      B[j][i]=B[i][j];
      ModTot=ModTot+B[i][j];
      if (A[i][j]>0.)
        ModMax=ModMax+B[i][j];
      A[j][i]=0.;  // La partie supérieure droite marquera à 1 les paires d'éléments réunis dans au moins une classe
    }
  }
}

int Louv1()
{
  int i, j, k, flag=1, fflag=0, NbTrans=0, OldC, NewC;
  float  VarMax;

  // Calcul de la contribution de chaque element a chaque classe
  for (i=0; i<N; i++)
  {
    for (k=1; k<=NbClas; k++)
      Var[i][k]=0.;
    for (j=0; j<N; j++)
    {
      k=Part[j];
      Var[i][k]=Var[i][k]+B[i][j];
    }
  }

  // fflag=1 si au moins 1 transfert dans cet appel de la procédure
  while(flag>0)
  {
    flag=0;  // transfert dans cette boucle
    for (i=0; i<N; i++)
    {
      OldC=Part[i];
      VarMax=Var[i][OldC];
      NewC=0;
      for (k=1; k<=NbClas; k++)
      {
        if (Var[i][k]>VarMax)
        {
          VarMax=Var[i][k];
          NewC=k;
        }
      }
      if (NewC>0 || VarMax<0.)
      {
        flag=1;
        fflag=1;
      }
      else
        continue;
      NbTrans++;
      for (j=0; j<N; j++)
        Var[j][OldC]=Var[j][OldC]-B[j][i];
      if (VarMax<0.)
      {
        NbClas++;
        NewC=NbClas;
        for (j=0; j<N; j++)
          Var[j][NewC]=0.;
      }
      Part[i]=NewC;
      for (j=0; j<N; j++)
        Var[j][NewC]=Var[j][NewC]+B[j][i];
    }
  }
  return fflag;
}

// Renumérote les classes et Calcule la modularité
void Renum()
{
  int i, j1, j2, k, kk=0, k1, k2, card;
  double  Mod, SumMod, MMod;

  for (k=0; k<=NbClas; k++)
    Kard[k]=0;
  for (i=0; i<N; i++)
    Kard[Part[i]]++;
  SumMod=0.; // Modularite du graphe non ponderé
  for (k=1; k<=NbClas; k++)
  {
    if (Kard[k]==0)
      continue;
    else
      kk++;
    card=0;
    Mod=0.;
    for (i=0; i<N; i++)
    {
      if (Part[i]==k)
      {
        Clas[card]=i;
        card++;
      }
    }
    for (i=0; i<card; i++)
      Part[Clas[i]]=kk;
    if (card<2)
      continue;
    for (k1=1; k1<card; k1++)
    {
      j1=Clas[k1];
      for (k2=0; k2<k1; k2++)
      {
        j2=Clas[k2];
        Mod=Mod+B[j1][j2];
      }
    }
    SumMod=SumMod+Mod;
  }
  NbClas=kk;
  MMod=SumMod/Na/Na/2-SumDg2/Na/Na/4;
}

int Louv2()
{
  int i, j, k,kk, flag=1,fflag=0, OldC, NewC;
  float  VarMax;

  // poids des connections entre classes
  for (k=0; k<=NbClas; k++)
  {
    for (kk=0; kk<=NbClas; kk++)
      A[k][kk]=0.;
  }
  for (i=1; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      A[Part[i]][Part[j]] = A[Part[i]][Part[j]]+B[i][j];
      A[Part[j]][Part[i]] = A[Part[j]][Part[i]]+B[i][j];
    }
  }

  for (k=1; k<=NbClas; k++)
  {
    Clas[k]=k;
    A[k][k]=0.;
    for (kk=1; kk<=NbClas; kk++)
      Var[k][kk]=A[k][kk];
  }

  // on fusionne des qu'il y a une connection > 0 entre classes
  while(flag>0)
  {
    flag=0;
    for (i=1; i<=NbClas; i++)
    {
      OldC=Clas[i];
      VarMax=Var[i][OldC];
      NewC=0;
      for (k=1; k<=NbClas; k++)
      {
        if (Var[i][k]>VarMax)
        {
          VarMax=Var[i][k];
          NewC=k;
        }
      }
      if (VarMax<0.)
      {
        NbClas++;
        NewC=NbClas;
        for (k=1; k<NbClas; k++)
        {
          Var[k][NbClas]=0.;
          Var[NbClas][k]=A[k][i];
        }
        Var[NbClas][OldC]=Var[i][OldC];
      }
      if (NewC>0)
      {
        flag=1;
        fflag=1;
      }
      else
        continue;
      Clas[i]=NewC; 
      // on deplace i de OldC a NewC ; mise a jour de Var
      for (j=1; j<=NbClas; j++) 
        Var[j][OldC]=Var[j][OldC]-A[j][i];
      for (j=1; j<=NbClas; j++)
        Var[j][NewC]=Var[j][NewC]+A[j][i];
    }
  }
  return fflag;
}

int TFit()
{
  int  i, NbPas=0, fflag=1;

  for (i=0; i<N; i++)
    Part[i]=i+1;
  NbClas=N;

  while (fflag)
  {
    NbPas++;
    fflag=Louv1();
    if (fflag==0)
      continue;
    Renum();
    // Y a-t-il des connections > 0 entre classes ?
    fflag=Louv2();
    if (fflag==0)
    {
      printf("\n");
      continue;
    }
    for (i=0; i<N; i++) Part[i]=Clas[Part[i]];
    Renum();
  }
  return NbPas;
}

int Transfert()
{
  int i,ii, j, k,kk, flag=1, NbTrans=0, OldC, NewC;
  float  VarMax, gain, gainmax;

  // Calcul de la contribution de chaque élément à chaque classe
  for (i=0; i<N; i++)
  {
    for (k=1; k<=NbClas; k++)
      Var[i][k]=0.;
    for (j=0; j<N; j++)
    {
      k=Part[j];
      Var[i][k]=Var[i][k]+B[i][j];
    }
    // meilleure affectation dans Clas[i]
    VarMax=-1.;
    for (k=1; k<=NbClas; k++)
    {
      if (Var[i][k]>VarMax)
      {
        VarMax = Var[i][k];
        kk=k;
      }
    }
    if (VarMax>0.)
      Clas[i]=kk;
    else
      Clas[i]=0;
    
    if (N<21)
    {
      printf("%*s : %2d, %2d : ",-NbCar,Et[i],Part[i],Clas[i]);
      for (k=1; k<=NbClas; k++)
        printf("%4.0f ",Var[i][k]);
      printf("\n");
    }
  }

  // fflag=1 si au moins 1 transfert dans cet appel de la proc?dure
  while (flag>0)
  {
    gainmax=0.;
    ii=-1;
    // cherche le meilleur transfert
    for (i=0; i<N; i++)
    {
      gain=Var[i][Clas[i]]-Var[i][Part[i]];
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
    OldC=Part[ii];
    NewC=Clas[ii];
    
    for (j=0; j<N; j++)
      Var[j][OldC] = Var[j][OldC] - B[j][ii];
    if (NewC==0)
    {
      NbClas++;
      NewC=NbClas;
    }
    Part[ii]=NewC;
    for (j=0; j<N; j++)
    {
      Var[j][NewC] = Var[j][NewC] + B[j][ii];
      if (Clas[j]==OldC) // on recalcule la meilleure affectation
      {
        VarMax=-1.;
        for (k=1; k<=NbClas; k++)
        {
          if (Var[j][k]>VarMax)
          {
            VarMax=Var[j][k];
            kk=k;
          }
        }
        if (VarMax>0.)
          Clas[j]=kk;
        else
          Clas[j]=0;
      }
      else if (Var[j][NewC]>Var[j][Clas[j]])
        Clas[j]=NewC;
    }
  }
  return NbTrans;
}

double Score()
{
  int i, j;
  double  Sc=0.;

  for (i=1; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      if (Part[i]==Part[j])
        Sc=Sc+B[i][j];
    }
  }
  return Sc;
}

void Around()
{
  int    i, j, ij, k=0, kmax, es=0, NbC, NbEs=N, NbTrans;
  double  ScIni, Sc;

  for (i=0; i<N; i++)
    Q[i]=Part[i];
  NbC=NbClas; // mémorise la meilleure partition
  ScIni=Score();
  printf("Sc Ini %.0f : ",ScIni);
  if (NbEs>500)
    NbEs=500;
  while (es<NbEs)
  {
    for (i=0; i<N; i++)
      Part[i]=Q[i]; // On repart de la meilleure
    NbClas=NbC;
    if (NbClas==1)
      return;
    k=0;
    kmax=1+(1.*rand()/RAND_MAX)*N/NbClas; // nb de swapping diminue avec les essais
    while (k<kmax)
    {
      i=(1.0*rand()/RAND_MAX)*N; // i indice au hasard entre 0 et N-1
      j=(1.0*rand()/RAND_MAX)*N; // j indice au hasard entre 0 et N-1
      if (Part[i]!=Part[j])
      {
        ij=Part[i];
        Part[i]=Part[j];
        Part[j]=ij;
        k++;
      }
    }

    NbTrans=Transfert();
    Sc=Score();
    if (Sc>ScIni)
    {
      for (i=0; i<N; i++)
        Q[i]=Part[i];
      ScIni=Sc;
      NbC=NbClas;
      printf("%d %.0f : ",kmax,Sc);
      es=0;
    }
    else
      es++;
  }
  printf("\n\n");
  NbClas=NbC;
  for (i=0; i<N; i++)
    Part[i]=Q[i];
  return;
}

void SaveClas(const char* FichS)
{
  int i, k, kk=0, Kmin=0;

  //printf("Minimum cardinality of selected classes (0=all) ");
  //scanf("%d",&Kmin); //TODO: could be a parameter
  Kmin = 0;  
  FILE* FichClas = fopen(FichS,"w"); assert(FichClas != NULL);
  
  for (k=0; k<=NbClas; k++)
    Kard[k]=0;
  for (i=0; i<N; i++)
    Kard[Part[i]]++;
  for (k=1; k<=NbClas; k++)
  {
    if (Kard[k]>=Kmin)
      kk++;
  }
  fprintf(FichClas," %d\n",kk);

  for (k=1; k<=NbClas; k++)
  {
    if (Kard[k]<Kmin)
      continue;
    fprintf(FichClas," %d\n",Kard[k]);
    for (i=0; i<N; i++)
    {
      if (Part[i]==k)
        fprintf(FichClas,"%s ",Et[i]);
    }
    fprintf(FichClas,"\n\n");
  }
  fprintf(FichClas,"\n\n");
  fclose(FichClas);
}

//X = chemin vers fichier... / out: fichier .clas
void tfit_core(const char* X, const char* out)
{
  int i, j,  NbPas=0, ClasMax;
  float  val;

  // calcule une partition qui optimise un critère de modularité
  // Algorithme de Transfert-Fusion itéré
  // Avec sauvegarde du fichier Clas

  LecGraph(X);

  printf("Nb. vertices %d, edges %d\n",N,Na);
  printf("percentage of edges in the graph %.4f\n\n",2.*Na/N/(N-1));
  
  A = malloc((N) * sizeof(float *));    // Graphe pondere
  B = malloc((N) * sizeof(float *));    // Matrice des modularites
  Var = malloc((N) * sizeof(float *));  // Variation en cas de fusion

  Clas = malloc((N+1) * sizeof(short));      // Classe de chaque element 
  Sum = malloc((N) * sizeof(float));      // Somme des poids des aretes en chaque sommet
  Kard = malloc((N+1) * sizeof(short));    // Nb. d'elements dans les classes courantes
  Part = malloc((N) * sizeof(short));      // Partition courante et finale 
  PartMax = malloc((N) * sizeof(short));      // Partition courante et finale 
  Q = malloc((N) * sizeof(short));      // subdivision d'une classe 

  assert (A != NULL  && B != NULL && Var != NULL);
  assert (Sum != NULL && Clas != NULL && Kard != NULL && Part != NULL && PartMax != NULL && Q != NULL);
  for (i=0; i<N; i++)
  {
    A[i] = malloc ( N * sizeof(float) );
    B[i] = malloc ( N * sizeof(float) );
    Var[i] = malloc ( (N+1) * sizeof(float) );
    assert (A[i] != NULL  && B[i] != NULL && Var[i] != NULL);
    Var[i][i]=0.;
    A[i][i]=0.;
  }
  // Le graphe donné ponderé par 1
  for (i=1; i<N; i++)
  {
    for (j=0; j<i; j++)
    {
      if (T[i][j]==1)
        val=1.;
      else
        val=0.;
      A[i][j]=val;
      A[j][i]=val;
    }
  }

  MatrixMod(1.0);
  NbPas=TFit();
  Around();
  printf("Classes after stochastic transfers\n\n");
  ClasOut();
  ClasMax=NbClas;

  SaveClas(out);
  printf("C'est fini\n");
  free(Dg);
}
