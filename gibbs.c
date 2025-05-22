#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_GENE_NUM 20 /*与えられるプロモータ領域の最大遺伝子数*/
#define MOTIF_LENGTH 7 //モチーフ長
#define CHARACTER_NUM 4 //文字の種類数

char character_type[CHARACTER_NUM] = {'A', 'C', 'G', 'T'}; // 塩基をchar型からint型に変換する配列
double g_base_fre_table[CHARACTER_NUM] = {7519429.0, 4637676.0, 4637676.0, 7519429.0}; //塩基出現頻度
double g_q[CHARACTER_NUM]={0.0}; //バックグラウンドの出現確率
double g_score[CHARACTER_NUM][MOTIF_LENGTH]={0.0}; //スコア行列
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

//文字の種類の判定
int check_char(char x){
  int type = 0;
  for(type = 0; type < CHARACTER_NUM; type++){
    if(x == character_type[type]){break;}
  }
  return type;
}

//出現位置の初期化
void init_pos(int matrix_n){
  srand((unsigned)time(NULL)); //乱数の初期化

  for(int gene_i = 0; gene_i < matrix_n; gene_i++){
    int pro_length = cal_matrix_length(g_pro[gene_i].seq); //プロモーター配列の長さ
    int rand_limit = pro_length - matrix_n + 1;
    g_pro[gene_i].pos = rand()%rand_limit;  //モチーフ出現位置のを乱数で決定
  }
}

//バックグラウンドの出現確率の計算
void cal_back_q(void){
  double sum_q = 0.0;
  for(int i = 0; i < CHARACTER_NUM; i++){
    sum_q += g_base_fre_table[i];
  }
  for(int i = 0; i < CHARACTER_NUM; i++){
    g_q[i] = g_base_fre_table[i] / sum_q;
    printf("%f ", g_q[i]);
    printf("\n");
  }
}

//モチーフ位置からのスコア行列の計算
void cal_motif_score(int selected_num, int matrix_n){

  //スコア行列の初期化
  for(int i = 0; i < CHARACTER_NUM; i++){
    for(int j = 0; j < MOTIF_LENGTH; j++){
      g_score[i][j] = 0;
    }
  }

  //出現頻度の計算
  for(int gene_i = 0; gene_i < matrix_n; gene_i++){
    if(gene_i == selected_num){continue;} //選ばれた1本の場合計算しない

    int start = g_pro[gene_i].pos; //モチーフの位置
    for(int i = start; i < g_pro[gene_i].pos + MOTIF_LENGTH; i++){
        g_score[check_char(g_pro[gene_i].seq[i])][i - start]++; //出現回数を加算  
    }
  }
  
  //スコアの表示
  for(int i = 0; i < CHARACTER_NUM; i++){
    printf("\n");
    for(int j = 0; j < MOTIF_LENGTH; j++){
      printf("%5.0f ", g_score[i][j]);
    }
  }

  //スコア行列の計算
  for(int i = 0; i < CHARACTER_NUM; i++){
    printf("\n");
    for(int j = 0; j < MOTIF_LENGTH; j++){
      g_score[i][j]++; //疑似頻度1を加算
      g_score[i][j] /= (matrix_n - 1); //出現確率の計算
      g_score[i][j] /= g_q[i]; //スコア行列の計算
      //printf("%1.5f ", g_score[i][j]);
    }
  }
}

//ギブスサンプリングによる結合部位の発見
int gibbs_scan(int matrix_n){
  init_pos(matrix_n);

  for(int gene_i = 0; gene_i < matrix_n; gene_i++){ //プロモータの中から一本選ぶ
    cal_motif_score(gene_i, matrix_n); //選んだ1本以外の配列からスコア行列を作成

  }


  
}

int main(int argc, char* argv[]){
  int gene_num = read_promoter(argv[1]);  //1番目の引数で指定した遺伝子のプロモータ領域を読み込む
  cal_back_q();
  gibbs_scan(gene_num);


  return 0;
}