#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_GENE_NUM 20 /*与えられるプロモータ領域の最大遺伝子数*/
#define MOTIF_LENGTH 7 //モチーフ長

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
  int pos;
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//プロモーター配列の読み込み
int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1);
  }

  while(fscanf(fp, "%s", buffer) != EOF){
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0';
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); 
    }else{
      strcpy(g_pro[gene_num].seq,buffer);
      gene_num++;
    }    
  }
  return gene_num;
}

//配列の長さの計算
int cal_matrix_length(char matrix[]){
  int length = 0;
  for(length = 0; length < BUFSIZE; length++){
        if(matrix[length]=='\0'){break;}
      }
  return length;
}

//出現位置の初期化
void init_pos(int n){
  srand((unsigned)time(NULL)); //乱数の初期化

  for(int gene_i = 0; gene_i < n; gene_i++){
    int pro_length = cal_matrix_length(g_pro[gene_i].seq); //プロモーター配列の長さ
    int rand_limit = pro_length - n + 1;
    g_pro[gene_i].pos = rand()%rand_limit;  //モチーフ出現位置のを乱数で決定
  }
}

int gibbs_scan(int n){
  init_pos(n);

  
}

int main(int argc, char* argv[]){
  int gene_num = read_promoter(argv[1]);  //1番目の引数で指定した遺伝子のプロモータ領域を読み込む

  gibbs_scan(gene_num);
  init_pos(gene_num);
  for(int i=0; i<gene_num; i++){
    printf("%d\n", g_pro[i].pos);
  }
  return 0;
}