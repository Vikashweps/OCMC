#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// Полином x^3 + x + 1  → 3-битный CRC
int poly[4] = {1, 0, 1, 1};  // x^3 + x^1 + x^0

// Функция вычисления CRC
void calculateCRC(int data[], int data_bits, int crc[]) {
    for(int i = 0; i < data_bits; i++) {
        int bit = data[i];
        int msb = crc[0];
        int fb = msb ^ bit;
        
        // Сдвигаем CRC регистр влево на 1 бит (3 бита → сдвигаем 0->1, 1->2)
        for(int j = 0; j < 2; j++) 
            crc[j] = crc[j + 1];
        crc[2] = 0;  // Младший бит заполняем нулем
        
        // Если обратная связь = 1, применяем полином (без старшего бита x^3)
        if(fb == 1) {
            for(int j = 0; j < 3; j++) 
                crc[j] ^= poly[j + 1];  // poly[1], poly[2], poly[3]
        }
    }
}

// Функция проверки CRC
int checkCRC(int data[], int total_bits) {
    int crc[3] = {0};  // 3-битный регистр
    
    for(int i = 0; i < total_bits; i++) {
        int bit = data[i];
        int msb = crc[0];
        int fb = msb ^ bit;
        
        for(int j = 0; j < 2; j++) 
            crc[j] = crc[j + 1];
        crc[2] = 0;
        
        if(fb == 1) {
            for(int j = 0; j < 3; j++) 
                crc[j] ^= poly[j + 1];
        }
    }
    
    for(int i = 0; i < 3; i++) {
        if(crc[i] != 0) return 0;
    }
    return 1;
}

int main() {
    srand(time(NULL));
    
    printf(" CRC проверка\n\n");
    
    printf("1. 32 бита:\n");
    
    int data32[32] = {1,0,1,1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,0,0,1,0,0,0,1,1,0,1,0,0};
    int crc32[3] = {0};
    
    calculateCRC(data32, 32, crc32);
    
    printf("Данные: ");
    for(int i = 0; i < 32; i++) printf("%d", data32[i]);
    
    printf("\nCRC: ");
    for(int i = 0; i < 3; i++) printf("%d", crc32[i]);
    
    int data_with_crc32[35];  // 32 + 3 = 35
    for(int i = 0; i < 32; i++) data_with_crc32[i] = data32[i];
    for(int i = 0; i < 3; i++) data_with_crc32[32 + i] = crc32[i];
    
    printf("\nПроверка: %s\n\n", checkCRC(data_with_crc32, 35) ? "Ошибок нет" : "Ошибка!");
    
    printf("2. 250 бит:\n");
    
    int data250[250];
    for(int i = 0; i < 250; i++) data250[i] = rand() % 2;
    
    printf("Данные: ");
    for(int i = 0; i < 250; i++) printf("%d", data250[i]);

    int crc250[3] = {0};
    calculateCRC(data250, 250, crc250);
    
    printf("\nCRC: ");
    for(int i = 0; i < 3; i++) printf("%d", crc250[i]);
    printf("\nДанные измененные: ");
    
    int data_with_crc250[253];  // 250 + 3 = 253
    for(int i = 0; i < 250; i++) data_with_crc250[i] = data250[i];
    for(int i = 0; i < 3; i++) data_with_crc250[250 + i] = crc250[i];
    
    int good = 0, bad = 0;
    for(int i = 0; i < 253; i++) {
        int test[253];
        for(int j = 0; j < 253; j++) test[j] = data_with_crc250[j];
        test[i] = 1 - test[i];
        
        if(checkCRC(test, 253)) {
            bad++;
        } else {
            good++;
        }
        printf("%d", test[i]);
    }
    
    printf("\nОбнаружено ошибок: %d/%d (%.1f%%)\n", good, 253, (float)good / 253 * 100);
    return 0;
}