#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdbool.h>

#define BUFSIZE 1024 //ファイルから読み込む一行の最大文字数
#define MAX_SEQ_NUM 30 //一つの転写因子に対して与えられる結合部位配列の最大数
#define MAX_GENE_NUM 8 /*与えられるプロモータ領域の最大遺伝子数*/
#define CHARACTER_NUM 4 //文字の種類数
#define THRESHOLD 4.9//閾値
#define RAND_LIMIT 100 //生成する乱数の最大値+1
#define RAND_PRO_NUM 10 //ランダム配列の数
#define ERROR_THRESHOLD 0.01 //ランダム配列の出現確率の誤差における閾値

char character_type[CHARACTER_NUM] = {'A', 'C', 'G', 'T'}; // 塩基をchar型からint型に変換する配列
char g_motif[MAX_SEQ_NUM][BUFSIZE]; //転写因子の結合部位配列を保存する配列
double g_fre_table[CHARACTER_NUM][BUFSIZE]={0.0}; //頻度表
double g_base_fre_table[CHARACTER_NUM] = {7519429.0, 4637676.0, 4637676.0, 7519429.0}; //塩基出現頻度
double g_q[CHARACTER_NUM]={0.0}; //バックグラウンドの出現確率
double g_log_odds_score[CHARACTER_NUM][BUFSIZE]={0.0}; //対数オッズスコア表
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

//頻度表の計算
void make_fre_table(int seq_num, int motif_length){
    for(int i = 0; i < seq_num; i++){
        for(int j = 0; j < motif_length; j++){
            g_fre_table[check_char(g_motif[i][j])][j]++;
        }
    }

    //頻度表の表示
    for(int i = 0; i < CHARACTER_NUM; i++){
      printf("%c | ", character_type[i]);
      for(int j = 0; j < motif_length; j++){
        printf("%5.0f ", g_fre_table[i][j]);
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
  double p[CHARACTER_NUM][motif_length];
  for(int j = 0; j < motif_length; j++){
    double sum_p = 0.0;
    for(int i = 0; i < CHARACTER_NUM; i++){
      sum_p += g_fre_table[i][j];
    }
    for(int i = 0; i < CHARACTER_NUM; i++){
      p[i][j] = g_fre_table[i][j] / sum_p;
    }
  }

  //バックグラウンドの出現確率の計算
  double sum_q = 0.0;
  for(int i = 0; i < CHARACTER_NUM; i++){
    sum_q += g_base_fre_table[i];
  }
  for(int i = 0; i < CHARACTER_NUM; i++){
    g_q[i] = g_base_fre_table[i] / sum_q;
    printf("%f ", g_q[i]);
    printf("\n");
  }

  //対数オッズスコアの計算
  for(int i = 0; i < CHARACTER_NUM; i++){
    for(int j = 0; j < motif_length; j++){
      g_log_odds_score[i][j] = log(p[i][j]/g_q[i]);
    }
  }

  //対数オッズスコアの表示
  for(int i = 0; i < CHARACTER_NUM; i++){
    printf("%c | ", character_type[i]);
    for(int j = 0; j < motif_length; j++){
        printf("%5.2f ", g_log_odds_score[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

//一致度の計算
double scan(int start, int motif_length, char pro[]){
  double score = 0.0;
  for(int i = start; i < start + motif_length; i++){
    score += g_log_odds_score[check_char(pro[i])][i - start];
  }
  return score;
}

//結合部位の探索
void hit(int motif_length, int gene_num){
  struct binding_site{
    double score;
    int pos;
  }b_site[gene_num];  //転写因子結合部位の候補のスコアと位置

  for(int gene_i = 0; gene_i < gene_num; gene_i++){

    //初期化
    b_site[gene_i].score = 0.0;
    int start = 0;
    
      while (g_pro[gene_i].seq[start + motif_length - 1] != '\0'){
      double score = scan(start, motif_length, g_pro[gene_i].seq);

      //出力
      if(THRESHOLD < score){
        b_site[gene_i].score = score;
        b_site[gene_i].pos = start;
        printf("pro:%s\n", g_pro[gene_i].name);
        printf("pos:%d\n", b_site[gene_i].pos + 1);
        printf("hit(");
        bool positive = true;
          for(int i = b_site[gene_i].pos; i < b_site[gene_i].pos + motif_length; i++){
            if(g_fre_table[check_char(g_pro[gene_i].seq[i])][i] == 1){positive = false;}
            printf("%c", g_pro[gene_i].seq[i]);
          }
        printf(")=%1.2f\n", b_site[gene_i].score);
        if(!positive){printf("偽陽性");}
        printf("\n");
      }
      start++;
    }
  }
}

//平均
double cal_ave(double score_data[][BUFSIZE], int row, int line){
  double sum = 0.0;
  double count = 0.0;
  for(int i = 0; i < row; i++){
    for(int j = 0; j < line; j++){
      sum += score_data[i][j];
      count++;
    }
  }
  return sum/count;
}

//標準偏差
double cal_sd(double score_data[][BUFSIZE], int row, int line, double ave){
  double sd = 0.0;
  double count = 0.0;
  for(int i = 0; i < row; i++){
    for(int j = 0; j < line; j++){
      sd += (score_data[i][j] - ave)*(score_data[i][j] - ave);
      count++;
    }
  }
  sd /= count;
  sd = sqrt(sd);
}

//乱数による文字の生成
char random_char(double r){
  double cumulative = 0.0;
    for(int i = 0; i < CHARACTER_NUM; ++i){
        cumulative += g_q[i];
        if(r < cumulative){
            return character_type[i];
        }
    }
  return character_type[CHARACTER_NUM-1]; //
}

//ランダム配列の出現確率の誤差計算
double cal_total_error(int count[], int length){
  double q_actual[CHARACTER_NUM] = {0.0};
  for(int i = 0; i < CHARACTER_NUM; i++){
    q_actual[i] = (double)count[i] / length;
  }

  double total_error = 0.0;
  for(int i = 0; i < CHARACTER_NUM; i++){
    total_error += fabs(q_actual[i] - g_q[i]);
  }
  return total_error;
}

//バックグラウンド確率に従ったランダム配列のスコア計算
void cal_random_score(int motif_length){
  srand((unsigned)time(NULL)); //乱数の初期化
  int start = 0;
  double rand_score[RAND_PRO_NUM][BUFSIZE]={0.0}; //ランダム配列のスコアを格納する配列
  int pro_length = cal_matrix_length(g_pro[0].seq); //作成するプロモータの長さ
  char rand_pro[pro_length+1]; //ランダムプロモータ配列を格納する配列

  int rand_num = 0; //生成したランダム配列の本数

  while(rand_num < RAND_PRO_NUM){
    int count_char[CHARACTER_NUM] = {0}; //ある文字の出現回数を計測する配列

    for(int i = 0; i < pro_length; i++){
      double tmp_rand = rand() / (RAND_MAX + 1.0); //[0,1)の乱数を生成
      char x = random_char(tmp_rand); 
      rand_pro[i] = x;

      for(int k = 0; k < CHARACTER_NUM; k++){
        if(x == character_type[k]){count_char[k]++;} //出現文字の回数を計測
      }
    }
    rand_pro[pro_length] = '\0'; //終端文字

    double error = cal_total_error(count_char, pro_length); 
    
    if(error < ERROR_THRESHOLD){ //ランダム配列が妥当なとき
      start = 0;
      while (rand_pro[start + motif_length - 1] != '\0'){ //全ての部分配列のスコア計算
        rand_score[rand_num][start] = scan(start, motif_length, rand_pro);
        start++;
      }

      rand_num++; //生成したランダム配列の個数を加算
    }
  }
    
  //平均値と標準偏差の計算
  double ave_score = cal_ave(rand_score, RAND_PRO_NUM, start);
  double sd_score = cal_sd(rand_score, RAND_PRO_NUM, start, ave_score);
  printf("ave:%f\nsd:%f\n", ave_score, sd_score);
}


int main(int argc, char* argv[]){
  int seq_num = read_multi_seq(argv[1]); //１番目の引数で指定した転写因子の複数の結合部位配列を読み込む
  int gene_num = read_promoter(argv[2]); //２番目の引数で指定した遺伝子のプロモータ領域を読み込む
  int motif_length = cal_matrix_length(g_motif[0]); //転写因子の結合部位領域配列の長さの計算

  make_fre_table(seq_num, motif_length); //頻度表の計算
  make_odds_score(motif_length); //対数オッズスコアの計算
  hit(motif_length, gene_num); //転写因子結合部位の予測
  cal_random_score(motif_length); //ランダム配列のスコアの計算
  
  return 0;
}