#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define CHARACTER_NUM 4 //文字の種類数
#define THRESHOLD 6.5 //閾値
#define RAND_LIMIT 100 //生成する乱数の最大値+1
#define RAND_PRO_NUM 10 //ランダム配列の数

enum dna {A,C,G,T}; //A, C, G, Tのchar型とint型の対応
char character_type[CHARACTER_NUM] = {'A', 'C', 'G', 'T'}; // 塩基をchar型からint型に変換する配列
char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列
int g_fre_table[CHARACTER_NUM][BUFSIZE]={0}; //頻度表
float g_base_fre_table[CHARACTER_NUM] = {7519429, 4637676, 4637676, 7519429}; //塩基出現頻度
float g_q[CHARACTER_NUM]; //バックグラウンドの出現確率
float g_log_odds_score[CHARACTER_NUM][BUFSIZE]; //対数オッズスコア表
struct promoter{
  char name[BUFSIZE];
  char seq[BUFSIZE];
}g_pro[MAX_GENE_NUM]; //遺伝子のプロモータ領域を保存する構造体

//転写因子の結合部位領域配列群ファイルの読み込み
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
  return seq_num; //結合部位領域配列の数
}

//プロモーター配列ファイルの読み込み
int read_promoter(char *filename){
  int gene_num = 0;  
  char buffer[BUFSIZE];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL){
    printf("scorefile open error.\n");
    exit(1); //ファイルが開けなかった場合プログラムを終了
  }

  while(fscanf(fp, "%s", buffer) != EOF){ //プログラムから一行ずつ読み込む
    if(buffer[strlen(buffer)-1]=='\n'){
      buffer[strlen(buffer)-1]='\0'; //改行を切り落とす
    }
    
    if(buffer[0]=='>'){
      strcpy(g_pro[gene_num].name,buffer+1); //ORF名を保存
    }else{
      strcpy(g_pro[gene_num].seq,buffer); //プロモーター配列を保存
      gene_num++; 
    }    
  }
  return gene_num; //遺伝子の数を返す
}

//配列の長さの計算
int matrix_length(void){
  int length = 0;
  for(length = 0; length < BUFSIZE; length++){
        if(g_motif[0][length]=='\0'){break;}
      }
  return length;
}

//プロモーター領域の長さ
int cal_pro_length(void){
  int length = 0;
  for(length = 0; length < BUFSIZE; length++){
    if(g_pro[0].seq[length]=='\0'){break;}
  }
  return length;
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

    //頻度表の表示
    for(int i = 0; i < CHARACTER_NUM; i++){
      printf("%c ", character_type[i]);
      for(int j = 0; j < motif_length; j++){
        printf("%3d ", g_fre_table[i][j]);
      }
      printf("\n");
    }
    printf("\n");
}

//対数オッズスコア行列の計算
void make_odds_score(int motif_length){
  //疑似頻度1を加える
  for(int i = 0; i < CHARACTER_NUM; i++){
    for(int j = 0; j < motif_length; j++){
         g_fre_table[i][j]++;
    }
  }

  //塩基の出現確率の計算
  float p[CHARACTER_NUM][motif_length];
  for(int j = 0; j < motif_length; j++){
    float sum_p = 0;
    for(int i = 0; i < CHARACTER_NUM; i++){
      sum_p += g_fre_table[i][j];
    }
    for(int i = 0; i < CHARACTER_NUM; i++){
      p[i][j] = g_fre_table[i][j] / sum_p;
    }
  }

  //バックグラウンドの出現確率の計算
  float sum_q = 0;
  for(int i = 0; i < CHARACTER_NUM; i++){
    sum_q += g_base_fre_table[i];
  }
  for(int i = 0; i < CHARACTER_NUM; i++){
    g_q[i] = g_base_fre_table[i] / sum_q;
  }

  //対数オッズスコアの計算
  for(int i = 0; i < CHARACTER_NUM; i++){
    for(int j = 0; j < motif_length; j++){
      g_log_odds_score[i][j] = log(p[i][j]/g_q[i]);
    }
  }

  //対数オッズスコアの表示
  for(int i = 0; i < CHARACTER_NUM; i++){
    for(int j = 0; j < motif_length; j++){
        printf("%2.2f ", g_log_odds_score[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

//一致度の計算
float scan(int start, int motif_length, char pro[]){
  float score = 0;
  for(int i = start; i < start + motif_length; i++){
          if(pro[i]=='A'){score += g_log_odds_score[A][i - start];}
          else if(pro[i]=='C'){score += g_log_odds_score[C][i - start];}
          else if(pro[i]=='G'){score += g_log_odds_score[G][i - start];}
          else if(pro[i]=='T'){score += g_log_odds_score[T][i - start];}
        }
  return score;
}

//結合部位の探索
void hit(int motif_length, int gene_num){
  struct binding_site{
    float score;
    int pos;
  }b_site[gene_num];  //転写因子結合部位の候補のスコアと位置
  for(int gene_i = 0; gene_i < gene_num; gene_i++){

    //初期化
    b_site[gene_i].score = 0;
    int start = 0;
    int b_num = 0;
    
      while (g_pro[gene_i].seq[start + motif_length - 1] != '\0'){
      float score = scan(start, motif_length, g_pro[gene_i].seq);

      //出力
      if(THRESHOLD < score){
        b_site[gene_i].score = score;
        b_site[gene_i].pos = start;
        printf("pro:%s\n", g_pro[gene_i].name);
        printf("pos:%d\n", b_site[gene_i].pos + 1);
        printf("hit(");
          for(int i = b_site[gene_i].pos; i < b_site[gene_i].pos + motif_length; i++){
            printf("%c",g_pro[gene_i].seq[i]);
          }
        printf(")=%1.2f\n", b_site[gene_i].score);
        printf("\n");
      }
      start++;
    }
  }
}

//平均
float cal_ave(float score_data[][BUFSIZE], int row, int line){
  float sum = 0;
  float count = 0;
  for(int i = 0; i < row; i++){
    for(int j = 0; j < line; j++){
      sum += score_data[i][j];
      count++;
    }
  }
  return sum/count;
}

//標準偏差
float cal_sd(float score_data[][BUFSIZE], int row, int line, float ave){
  float sd = 0;
  int count = 0;
  for(int i = 0; i < row; i++){
    for(int j = 0; j < line; j++){
      sd += (score_data[i][j] - ave)*(score_data[i][j] - ave);
      count++;
    }
  }
  sd /= count;
  sd = sqrt(sd);
}

//バックグラウンド確率に従ったランダム配列のスコア計算
void cal_random_score(int motif_length){
  srand((unsigned)time(NULL));
  int start = 0;
  float rand_score[RAND_PRO_NUM][BUFSIZE];
  int pro_length = cal_pro_length();
  char rand_pro[pro_length];
  for(int rand_num = 0; rand_num < RAND_PRO_NUM; rand_num++){
    for(int i = 0; i < pro_length; i++){
      int tmp_rand = rand()%RAND_LIMIT;
      if(tmp_rand<=30){
        rand_pro[i] = 'A';
      }
      else if(tmp_rand>=31&&tmp_rand<=49){
        rand_pro[i] = 'C';
      }
      else if(tmp_rand>=50&&tmp_rand<=68){
        rand_pro[i] = 'G';
      }
      else if(tmp_rand>=69){
        rand_pro[i] = 'T';
      }
    }
    start = 0;
    while (rand_pro[start + motif_length - 1] != '\0'){
      rand_score[rand_num][start] = scan(start, motif_length, rand_pro);
      //printf("%f \n",rand_score[rand_num][start]);
      start++;
    } 
  }

  float ave_score = cal_ave(rand_score, RAND_PRO_NUM, start);
  float sd_score = cal_sd(rand_score, RAND_PRO_NUM, start, ave_score);
  printf("ave:%f\nsd:%f\n", ave_score, sd_score);
}


int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  int gene_num = read_promoter(argv[2]);  //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  int motif_length = matrix_length();　//転写因子の結合部位領域配列の長さの計算

  make_fre_table(seq_num, motif_length); //頻度表の計算
  make_odds_score(motif_length); //対数オッズスコアの計算
  hit(motif_length, gene_num); //転写因子結合部位の予測
  cal_random_score(motif_length); //ランダム配列のスコアの計算

  return 0;
}