#include <stdio.h>

#define LENGTH 31
#define N 5

// Генерация M-последовательностей
void generate_x_sequence(int reg[], int output[]) {
    for (int i = 0; i < LENGTH; i++) {
        output[i] = reg[N - 1];
        int new_bit = reg[3] ^ reg[4];
        for (int j = N - 1; j > 0; j--) {
            reg[j] = reg[j - 1];
        }
        reg[0] = new_bit;
    }
}

void generate_y_sequence(int reg[], int output[]) {
    for (int i = 0; i < LENGTH; i++) {
        output[i] = reg[N - 1];
        int new_bit = reg[2] ^ reg[4];
        for (int j = N - 1; j > 0; j--) {
            reg[j] = reg[j - 1];
        }
        reg[0] = new_bit;
    }
}

// Циклический сдвиг
void cyclic_shift(int original[], int shifted[], int shift) {
    for (int i = 0; i < LENGTH; i++) {
        shifted[i] = original[(i + shift) % LENGTH];
    }
}

// Вычисление корреляции
double correlation(int seq1[], int seq2[]) {
    int sum = 0;
    for (int i = 0; i < LENGTH; i++) {
        sum += (seq1[i] == seq2[i]) ? 1 : -1;
    }
    return (double)sum / LENGTH;
}

int main() {
    // Первая последовательность Голда (x=12, y=19)
    int x_reg1[N] = {0, 1, 1, 0, 0};
    int y_reg1[N] = {1, 0, 0, 1, 1};
    int x_seq1[LENGTH], y_seq1[LENGTH], gold_seq1[LENGTH];

    generate_x_sequence(x_reg1, x_seq1);
    generate_y_sequence(y_reg1, y_seq1);
    
    for (int i = 0; i < LENGTH; i++) {
        gold_seq1[i] = x_seq1[i] ^ y_seq1[i];
    }

    printf("Первая последовательность Голда (x=12, y=19):\n");
    for (int i = 0; i < LENGTH; i++) {
        printf("%d", gold_seq1[i]);
    }
    printf("\n\n");

    // Таблица автокорреляции для первой последовательности
    printf("Таблица автокорреляции первой последовательности:\n");
    printf("Сдвиг | бит 1 | бит 2 | бит 3 | Автокорреляция\n");
    printf("------|-------|-------|-------|----------------\n");

    for (int shift = 0; shift < LENGTH; shift++) {
        int shifted_seq[LENGTH];
        cyclic_shift(gold_seq1, shifted_seq, shift);
        double R_tau = correlation(gold_seq1, shifted_seq);
        
        printf("%5d | %5d  | %5d  | %5d  | %+8.3f\n", 
               shift, shifted_seq[0], shifted_seq[1], shifted_seq[2], R_tau);
    }
    printf("\n");

    // Вторая последовательность Голда (x=13, y=14)
    int x_reg2[N] = {0, 1, 1, 0, 1};
    int y_reg2[N] = {0, 1, 1, 1, 0};
    int x_seq2[LENGTH], y_seq2[LENGTH], gold_seq2[LENGTH];

    generate_x_sequence(x_reg2, x_seq2);
    generate_y_sequence(y_reg2, y_seq2);
    
    for (int i = 0; i < LENGTH; i++) {
        gold_seq2[i] = x_seq2[i] ^ y_seq2[i];
    }

    printf("Вторая последовательность Голда (x=13, y=14):\n");
    for (int i = 0; i < LENGTH; i++) {
        printf("%d", gold_seq2[i]);
    }
    printf("\n\n");

    // Взаимная корреляция
    double cross_correlation = correlation(gold_seq1, gold_seq2);
    printf("Взаимная корреляция между последовательностями: %+8.3f\n", cross_correlation);

    return 0;
}