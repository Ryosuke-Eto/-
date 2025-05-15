#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define BASE 4

enum dna {A,C,G,T};
char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列
int g_fre_table[BASE][BUFSIZE]={0}; //頻度表
float base_fre_table[BASE] = {7519429, 4637676, 4637676, 7519429}; //塩基出現頻度
float log_odds_score[BASE][BUFSIZE]; //対数オッズスコア表

struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//グローバル変数はローカル変数と区別するため、全部大文字にするかg_を先頭につけるのが一般的

int read_multi_seq(char* filename){
  int seq_num = 0;
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("motif_region_file open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    strcpy(g_motif[seq_num],buffer); //結合部位配列を保存
    seq_num++;
  }
  return seq_num;
}

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

//頻度表の計算
void make_fre_table(int seq_num, int motif_length){
    for(int i = 0; i < seq_num; i++){
        for(int j = 0; j < motif_length; j++){
            if(g_motif[i][j]=='A'){g_fre_table[A][j]++;}
            else if(g_motif[i][j]=='C'){g_fre_table[C][j]++;}
            else if(g_motif[i][j]=='G'){g_fre_table[G][j]++;}
            else if(g_motif[i][j]=='T'){g_fre_table[T][j]++;}
        }
    }
}

//対数オッズスコア行列の計算
void make_odds_score(int motif_length){
  //疑似頻度1を加える
  for(int i = 0; i < BASE; i++){
    for(int j = 0; j < motif_length; j++){
         g_fre_table[i][j]++;
    }
  }

  //塩基の出現確率の計算
  float p[BASE][motif_length];
  float sum_p = 0;
  for(int j = 0; j < motif_length; j++){
    for(int i = 0; i < BASE; i++){
      sum_p += g_fre_table[i][j];
    }
    for(int i = 0; i < BASE; i++){
      p[i][j] = g_fre_table[i][j] / sum_p;
    }
  }

  //バックグラウンドの出現確率の計算
  float q[BASE];
  float sum_q = 0;
  for(int i = 0; i < BASE; i++){
    sum_q += base_fre_table[i];
  }
  for(int i = 0; i < BASE; i++){
    q[i] = base_fre_table[i] / sum_q;
  }

  //対数オッズスコアの計算
  for(int i = 0; i < BASE; i++){
    for(int j = 0; j < motif_length; j++){
      log_odds_score[i][j] = log(p[i][j]/q[i]);
    }
  }
}

  //結合部位の探索
  void hit(int motif_length, int gene_num){
    struct score{
      float score;
      int pos;
    }match[gene_num];
    for(int gene_i = 0; gene_i < gene_num; gene_i++){
      
      //プロモーターの長さ
      int pro_length = 0;
      for(pro_length = 0; pro_length < BUFSIZE; pro_length++){
        if(g_pro[gene_i].seq[pro_length]=='\0'){break;}
      }

      //ゲノム配列上の結合部位の探索
      match[gene_i].score = 0;
      int end = pro_length - motif_length;
      for(int start = 0; start < end; start++){
        int score = 0;
        for(int i = start; i < start + motif_length; i++){
          score += log_odds_score[g_pro[gene_i].seq[i]][i-start];
        }

        if(match[gene_i].score < score){
          match[gene_i].score = score;
          match[gene_i].pos = start;
        }
      }

      printf("pro:%s\n", g_pro[gene_i].name);
      printf("pos:%d\n", match[gene_i].pos);
      printf("hit(");
      for(int i = match[gene_i].pos; i < match[gene_i].pos + motif_length; i++){
        printf("%c", g_pro[gene_i].seq[i]);
      }
      printf(")= %f\n", match[gene_i].score);
      printf("\n");
    }
  }

int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む

  //配列の長さの計算
  int motif_length = 0;
      for(motif_length = 0; motif_length < BUFSIZE; motif_length++){
        if(g_motif[0][motif_length]=='\0'){break;}
      }

  make_fre_table(seq_num, motif_length);

  //頻度表の表示
  for(int i = 0; i < BASE; i++){
    for(int j = 0; j < motif_length; j++){
        printf("%d ", g_fre_table[i][j]);
    }
    printf("\n");
  }
    
  make_odds_score(motif_length);
   
  //対数オッズスコアの表示
  for(int i = 0; i < BASE; i++){
    for(int j = 0; j < motif_length; j++){
        printf("%f ", log_odds_score[i][j]);
    }
    printf("\n");
  }

  hit(motif_length, gene_num);
  return 0;
}