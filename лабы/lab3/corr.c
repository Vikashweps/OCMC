#include <stdio.h>
#include <math.h>

int main() {
    double a[] = {6, 2, 8, -2, -4, -4, 1, 3};
    double b[] = {3, 6, 7, 0, -5, -4, 2, 5};
    double c[] = {-6, -1, -3, -9, 2, -8, 4, 1};
    
    // Обычная корреляция
    double ab = 0, ac = 0, bc = 0;
    double aa = 0, bb = 0, cc = 0;
    
    for(int i = 0; i < 8; i++) {
        ab += a[i] * b[i];
        ac += a[i] * c[i];
        bc += b[i] * c[i];
        aa += a[i] * a[i];
        bb += b[i] * b[i];
        cc += c[i] * c[i];
    }
    
    // Нормализованная корреляция
    double nab = ab / sqrt(aa * bb);
    double nac = ac / sqrt(aa * cc);
    double nbc = bc / sqrt(bb * cc);
    
    printf("\nCorrelation:\n");
    printf("  |   a  |   b  |   c\n");
    printf("-----------------------\n");
    printf("a | %4.0f | %4.0f | %4.0f\n", aa, ab, ac);
    printf("b | %4.0f | %4.0f | %4.0f\n", ab, bb, bc);
    printf("c | %4.0f | %4.0f | %4.0f\n\n", ac, bc, cc);
    
    printf("normal correlation:\n");
    printf("  |   a   |   b   |   c\n");
    printf("-------------------------\n");
    printf("a | %5.2f | %5.2f | %5.2f\n", 1.0, nab, nac);
    printf("b | %5.2f | %5.2f | %5.2f\n", nab, 1.0, nbc);
    printf("c | %5.2f | %5.2f | %5.2f\n", nac, nbc, 1.0);
    
    return 0;
}