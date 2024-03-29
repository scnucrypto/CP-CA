#include "wbspeck.h"

void Key_expand32(u16 key[4], u16 roundkey[22])
{
    u16 i, b = key[3];
    u16 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[i + 1];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R32(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}
void Key_expand64(u32 key[4], u32 roundkey[22])
{
    u32 i, b = key[3];
    u32 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[2 - i];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R64(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}
void Key_expand128(u64 key[4], u64 roundkey[22])
{
    u64 i, b = key[3];
    u64 a[SPECK_KEY_LEN - 1];

    for(i = 0; i < (SPECK_KEY_LEN - 1); i++)
    {
        a[i] = key[2 - i];
    }
    roundkey[0] = b;  
    for(i = 0; i < SPECK_ROUNDS - 1; i++)
    {
        R128(a[i % (SPECK_KEY_LEN - 1)], b, i);
        roundkey[i + 1] = b;
    }
}

u16 ROR32(u16 x, int r)
{
    return (x >> r) | (x << ((sizeof(u16) * 8) - r));
}
u16 ROL32(u16 x, int r)
{
    return (x << r) | (x >> ((sizeof(u16) * 8) - r));
}
u32 ROR64(u32 x, int r)
{
    return (x >> r) | (x << ((sizeof(u32) * 8) - r));
}
u32 ROL64(u32 x, int r)
{
    return (x << r) | (x >> ((sizeof(u32) * 8) - r));
}
u64 ROR128(u64 x, int r)
{
    return (x >> r) | (x << ((sizeof(u64) * 8) - r));
}
u64 ROL128(u64 x, int r)
{
    return (x << r) | (x >> ((sizeof(u64) * 8) - r));
}

void sub_DCA32_KeyCheck(u16 *key_candidate, int key_candidate_count, Aff32 Aff, int q, u16 *reduced_key, int *reduced_key_count)
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u32 state;
    
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key16, roundkey);

    int key_count = 0;

    u32 map[apt];
    int f = 0, g = 0; 
    FILE *fp = fopen("Result_SPECK32.txt", "a");
    printf("------Partial Key Verification-----\n");  
    fprintf(fp, "------Partial Key Verification-----\n"); 

    u16 key_guess[100] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < apt; x++) // adaptive inputs
        {
            x0 = (u16)(x) << ((q - 1) * 4);
            x1 = 0;
            
            y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
            y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
            
            z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
            z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

            state = (z0 << 16) | z1;
            map[x] = state;
            map[x] = affineU32(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[apt] = {0.0};
        int L_count[apt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < apt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb32; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < apt; x++)
                { 
                    if(map[x] & idM32[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < apt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count >= 28) //(l_count > k_count)
            {
                if(l_count > k_count) k_count = l_count;
                // knum = 0;
                key_guess[knum] = k;
                knum++;
                // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fp = fopen("Result_SPECK32.txt", "a");
                // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fclose(fp);
            }
            /*
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fp = fopen("Result_SPECK32.txt", "a");
                // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fclose(fp);
            }
            */
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < apt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "Verified key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        printf("Verified key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        reduced_key[k] = key_guess[k];
        key_count++;
    }
    *reduced_key_count = key_count;
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "Verified key count: %d\n", key_count);
    printf("Verified key count: %d\n", key_count);
    fclose(fp); 
}

void DCA32_KeyCheck(u16 *key_candidate, int key_candidate_count, Aff32 Aff)
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u32 state;
    
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key16, roundkey);

    int key_count = 0;

    int f = 0, g = 0;
    u32 map[pt]; 
    FILE *fp = fopen("Result_SPECK32.txt", "a");
    printf("------Final Key Verification-----\n");  
    fprintf(fp, "------Final Key Verification-----\n"); 

    u16 key_guess[100] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = (u16)(x) << 5;
            x1 = 1;
            
            y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
            y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
            
            z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
            z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

            state = (z0 << 16) | z1;
            map[x] = state;
            map[x] = affineU32(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int L_count[pt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb32; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(map[x] & idM32[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "simulation WB-SPECK32 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        printf("simulation WB-SPECK32 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        key_count++;
    }
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    fclose(fp); 
}
void DCA32()
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, i, j, s, ts, l, b, bit, q, len;
    u32 state;
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key16, roundkey);
    ////

    int f = 0, g = 0;
    Aff32 Aff, Aff_inv;
    genaffinepairM32(&Aff, &Aff_inv);
    u32 map[pt]; 
    
    FILE *fp; 
    u16 key_candidate[100] = {0};
    int key_candidate_count = 1;
    u16 temp_key_candidate[100] = {0};
    int temp_key_candidate_count = 0;
    
    for(q = 0; q < 4; q++)
    {
        memset(temp_key_candidate, 0, sizeof(temp_key_candidate));
        temp_key_candidate_count = 0;
        for(s = 0; s < key_candidate_count; s++)
        {
            u16 key_g = key_candidate[s];
            u16 key_guess[20] = {0};
            int k_count = 0;
            double k_max = 0.0;
            int knum = 0;
            for(kkey = 0; kkey <= kte; kkey++) // key, ktb, kte
            {
                k = (kkey << 4 * q) ^ key_g;
                for(x = 0; x < pt; x++) // adaptive inputs
                {
                    x0 = (u16)(x) << 4 * q;
                    x1 = 0;
                    
                    y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
                    y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
                    
                    z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
                    z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

                    state = (z0 << 16) | z1;
                    map[x] = state;
                    map[x] = affineU32(Aff, state); // affine encoding
                }
                double k_score = 0.0;
                double L_score[pt] = {0.0};
                int L_count[pt] = {0};
                int l_count = 0;
                double l_max = 0.0;
                for(l = 1; l < pt; l++)
                {  
                    double score = 0.0;
                    double j_max = 0.0;
                    int j_count = 0;
                    for(j = 0; j < tb32; j++) // samples of traces
                    {
                        int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                        for(x = 0; x < pt; x++)
                        { 
                            if(map[x] & idM32[j]) ts = 1; // traces
                            else ts = 0;
                            ////// correlation
                            f = ts;
                            if(f) Nf1++;
                            else Nf0++;
                            
                            g = xor[l & x];
                            if(g) Ng1++;
                            else Ng0++;
                            
                            if(f == 1 && g == 1) N11++;
                            else if(f == 0 & g == 0) N00++;
                            else if(f == 1 & g == 0) N10++;
                            else N01++;
                        }
                        if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                        else score = 0.0;
                        if(score > j_max)
                        {
                            j_count = 0;
                            j_max = score;
                            j_count++;
                        }
                        else if(fabs(score - j_max) < EPS)
                        {
                            j_count++;
                        }
                    }
                    L_score[l] = j_max;
                    L_count[l] = j_count;
                    if(L_score[l] > l_max) l_max = L_score[l];
                }
                k_score = l_max;
                if(fabs(k_score - k_max) < EPS)
                {
                    for(l = 1; l < pt; l++)
                    {
                        if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                        {
                            l_count += L_count[l];
                        }
                    }
                    if(l_count >= 28) //(l_count > k_count)
                    {
                        if(l_count > k_count) k_count = l_count;
                        // knum = 0;
                        key_guess[knum] = k;
                        knum++;
                        // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fp = fopen("Result_SPECK32.txt", "a");
                        // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fclose(fp);
                    }
                    /*
                    else if(l_count == k_count)
                    {
                        key_guess[knum] = k;
                        knum++;
                        // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fp = fopen("Result_SPECK32.txt", "a");
                        // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fclose(fp);
                    }
                    */
                }
                else if(k_score > k_max) 
                {
                    k_max = k_score; 
                    knum = 0;
                    for(l = 1; l < pt; l++)
                    {
                        if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                        {
                            l_count += L_count[l];
                        }
                    }
                    k_count = l_count;
                    key_guess[knum] = k;
                    knum++;
                }
            }
            fp = fopen("Result_SPECK32.txt", "a");
            printf("------Partial Key Guess-----\n");  
            fprintf(fp, "------Partial Key Guess-----\n");
            printf("------\n");  
            fprintf(fp, "------\n");  
            for(k = 0; k < knum; k++)
            {
                fprintf(fp, "No.%d nibble: recoverd key = %.2x, encoding count = %d, correlation: %f\n", q, key_guess[k], k_count, k_max);
                printf("No.%d nibble: recoverd key = %.2x, encoding count = %d, correlation: %f\n", q, key_guess[k], k_count, k_max);
                temp_key_candidate[temp_key_candidate_count] = key_guess[k];
                temp_key_candidate_count++;
            }
            printf("------\n"); 
            fprintf(fp, "------\n");  
            fclose(fp); 
        }
        len = sizeof(temp_key_candidate) / sizeof(temp_key_candidate[0]);
        memcpy(key_candidate, temp_key_candidate, len * sizeof(u16));
        key_candidate_count = temp_key_candidate_count;
        if(q > 0) sub_DCA32_KeyCheck(key_candidate, key_candidate_count, Aff, q, key_candidate, &key_candidate_count);    
    }
    
    fp = fopen("Result_SPECK32.txt", "a");
    printf("------Possible Key-----\n");  
    fprintf(fp, "------Possible Key-----\n");  
    for(i = 0; i < key_candidate_count; i++)
    {
        printf("%.2x\n", key_candidate[i]);
        fprintf(fp, "%.2x\n", key_candidate[i]);
    }
    printf("Count: %d\n\n", key_candidate_count);  
    fprintf(fp, "Count: %d\n\n", key_candidate_count); 
    fclose(fp); 

    DCA32_KeyCheck(key_candidate, key_candidate_count, Aff);
    
}

void sub_DCA64_KeyCheck(u32 *key_candidate, int key_candidate_count, Aff64 Aff, int q, u32 *reduced_key, int *reduced_key_count)
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u64 state;
    
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key32, roundkey);

    int key_count = 0;

    u64 map[apt]; 
    int f = 0, g = 0; 
    FILE *fp = fopen("Result_SPECK64.txt", "a");
    printf("------Partial Key Verification-----\n");  
    fprintf(fp, "------Partial Key Verification-----\n"); 

    u32 key_guess[100] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < apt; x++) // adaptive inputs
        {
            x0 = (u32)(x) << ((q - 1) * 4);
            x1 = 0;
            
            y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
            y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state = ((u64)(z0) << 32) | z1;
            map[x] = state;
            map[x] = affineU64(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[apt] = {0.0};
        int L_count[apt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < apt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb64; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < apt; x++)
                { 
                    if(map[x] & idM64[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < apt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count >= 56) //(l_count > k_count)
            {
                if(l_count > k_count) k_count = l_count;
                // knum = 0;
                key_guess[knum] = k;
                knum++;
                // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fp = fopen("Result_SPECK32.txt", "a");
                // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fclose(fp);
            }
            /*
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fp = fopen("Result_SPECK32.txt", "a");
                // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fclose(fp);
            }
            */
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < apt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "Verified key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        printf("Verified key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        reduced_key[k] = key_guess[k];
        key_count++;
    }
    *reduced_key_count = key_count;
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "Verified key count: %d\n", key_count);
    printf("Verified key count: %d\n", key_count);
    fclose(fp); 
}
void DCA64_KeyCheck(u32 *key_candidate, int key_candidate_count, Aff64 Aff)
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u64 state;
    
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key32, roundkey);

    int key_count = 0;

    int f = 0, g = 0;
    u64 map[pt]; 
    FILE *fp = fopen("Result_SPECK64.txt", "a");
    printf("------Final Key Verification-----\n");  
    fprintf(fp, "------Final Key Verification-----\n"); 

    u32 key_guess[100] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = (u32)(x) << 5;
            x1 = 1;
            
            y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
            y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state = ((u64)(z0) << 32) | z1;
            map[x] = state;
            map[x] = affineU64(Aff, state); // affine encoding
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int L_count[pt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb64; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(map[x] & idM64[j]) ts = 1; // traces
                    else ts = 0;
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "simulation WB-SPECK64 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        printf("simulation WB-SPECK64 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_count);
        key_count++;
    }
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "simulation WB-SPECK64 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK64 recoverd key count: %d\n", key_count);
    fclose(fp); 
}
void DCA64()
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, i, j, s, ts, l, b, bit, q, len;
    u64 state;
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key32, roundkey);
    ////

    int f = 0, g = 0;
    Aff64 Aff, Aff_inv;
    genaffinepairM64(&Aff, &Aff_inv);
    u64 map[pt]; 
    
    FILE *fp;

    u32 key_candidate[100] = {0};
    int key_candidate_count = 1;
    u32 temp_key_candidate[100] = {0};
    int temp_key_candidate_count = 0;
    for(q = 0; q < 8; q++)
    {
        memset(temp_key_candidate, 0, sizeof(temp_key_candidate));
        temp_key_candidate_count = 0;
        for(s = 0; s < key_candidate_count; s++)
        {
            u32 key_g = key_candidate[s];
            u32 key_guess[20] = {0};
            int k_count = 0;
            double k_max = 0.0;
            int knum = 0;
            for(kkey = ktb; kkey <= kte; kkey++) // key
            {
                k = (kkey << 4 * q) ^ key_g;
                for(x = 0; x < pt; x++) // adaptive inputs
                {
                    x0 = (u32)(x) << 4 * q;
                    x1 = 0;
                    
                    y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
                    y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
                    
                    z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
                    z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

                    state = ((u64)(z0) << 32) | z1;
                    map[x] = state;
                    map[x] = affineU64(Aff, state); // affine encoding
                }
                double k_score = 0.0;
                double L_score[pt] = {0.0};
                int L_count[pt] = {0};
                int l_count = 0;
                double l_max = 0.0;
                for(l = 1; l < pt; l++)
                {  
                    double score = 0.0;
                    double j_max = 0.0;
                    int j_count = 0;
                    for(j = 0; j < tb64; j++) // samples of traces
                    {
                        int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                        for(x = 0; x < pt; x++)
                        { 
                            if(map[x] & idM64[j]) ts = 1; // traces
                            else ts = 0;
                            ////// correlation
                            f = ts;
                            if(f) Nf1++;
                            else Nf0++;
                            
                            g = xor[l & x];
                            if(g) Ng1++;
                            else Ng0++;
                            
                            if(f == 1 && g == 1) N11++;
                            else if(f == 0 & g == 0) N00++;
                            else if(f == 1 & g == 0) N10++;
                            else N01++;
                        }
                        if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                        else score = 0.0;
                        if(score > j_max)
                        {
                            j_count = 0;
                            j_max = score;
                            j_count++;
                        }
                        else if(fabs(score - j_max) < EPS)
                        {
                            j_count++;
                        }
                    }
                    L_score[l] = j_max;
                    L_count[l] = j_count;
                    if(L_score[l] > l_max) l_max = L_score[l];
                }
                k_score = l_max;
                if(fabs(k_score - k_max) < EPS)
                {
                    for(l = 1; l < pt; l++)
                    {
                        if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                        {
                            l_count += L_count[l];
                        }
                    }
                    if(l_count >= 56) //(l_count > k_count)
                    {
                        if(l_count > k_count) k_count = l_count;
                        // knum = 0;
                        key_guess[knum] = k;
                        knum++;
                        // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fp = fopen("Result_SPECK32.txt", "a");
                        // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fclose(fp);
                    }
                    /*
                    else if(l_count == k_count)
                    {
                        key_guess[knum] = k;
                        knum++;
                        // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fp = fopen("Result_SPECK32.txt", "a");
                        // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fclose(fp);
                    }
                    */
                }
                else if(k_score > k_max) 
                {
                    k_max = k_score; 
                    knum = 0;
                    for(l = 1; l < pt; l++)
                    {
                        if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                        {
                            l_count += L_count[l];
                        }
                    }
                    k_count = l_count;
                    key_guess[knum] = k;
                    knum++;
                }
            }
            fp = fopen("Result_SPECK64.txt", "a");
            printf("------Partial Key Guess-----\n");  
            fprintf(fp, "------Partial Key Guess-----\n");
            printf("------\n");  
            fprintf(fp, "------\n");  
            for(k = 0; k < knum; k++)
            {
                fprintf(fp, "No.%d nibble: recoverd key = %.2x, encoding count = %d, correlation: %f\n", q, key_guess[k], k_count, k_max);
                printf("No.%d nibble: recoverd key = %.2x, encoding count = %d, correlation: %f\n", q, key_guess[k], k_count, k_max);
                temp_key_candidate[temp_key_candidate_count] = key_guess[k];
                temp_key_candidate_count++;
            }
            printf("------\n"); 
            fprintf(fp, "------\n"); 
            fclose(fp);  
        }
        len = sizeof(temp_key_candidate) / sizeof(temp_key_candidate[0]);
        memcpy(key_candidate, temp_key_candidate, len * sizeof(u32));
        key_candidate_count = temp_key_candidate_count;
        if(q > 0) sub_DCA64_KeyCheck(key_candidate, key_candidate_count, Aff, q, key_candidate, &key_candidate_count);
    }
    fp = fopen("Result_SPECK64.txt", "a");
    printf("------Possible Key-----\n");  
    fprintf(fp, "------Possible Key-----\n");  
    for(i = 0; i < key_candidate_count; i++)
    {
        printf("%.2x\n", key_candidate[i]);
        fprintf(fp, "%.2x\n", key_candidate[i]);
    }
    printf("Count: %d\n\n", key_candidate_count);  
    fprintf(fp, "Count: %d\n\n", key_candidate_count); 
    fclose(fp); 

    DCA64_KeyCheck(key_candidate, key_candidate_count, Aff);
}

void sub_DCA128_KeyCheck(u64 *key_candidate, int key_candidate_count, Aff128 Aff, int q, u64 *reduced_key, int *reduced_key_count)
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, j, s, ts, l, b, bit, i;
    u64 k, kkey;
    V128 state1;
    V128 state2;
    
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key64, roundkey);

    int key_count = 0;

    u64 map[apt][2]; 
    int f = 0, g = 0; 
    FILE *fp = fopen("Result_SPECK128.txt", "a");
    printf("------Partial Key Verification-----\n");  
    fprintf(fp, "------Partial Key Verification-----\n"); 

    u64 key_guess[200] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < apt; x++) // adaptive inputs
        {
            x0 = (u64)(x) << ((q - 1) * 4);
            x1 = 0;
            
            y1 = ROR128(ROL128(x0 - x1, 8) ^ x1, 3);
            y0 = ROL128((ROL128(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL128(y1, 3) ^ (ROR128(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR128((ROR128(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state1.V[0] = z0; 
            state1.V[1] = z1;
            MatMulVecM128(Aff.Mat, state1, &state2);
            map[x][0] = state2.V[0] ^ Aff.Vec.V[0];//add
            map[x][1] = state2.V[1] ^ Aff.Vec.V[1];
        }
        double k_score = 0.0;
        double L_score[apt] = {0.0};
        int L_count[apt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < apt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb128; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < apt; x++)
                { 
                    if(j < 64)
                    {
                        if(map[x][0] & idM64[j]) ts = 1; // traces
                        else ts = 0;
                    }
                    else
                    {
                        if(map[x][1] & idM64[j - 64]) ts = 1; // traces
                        else ts = 0;
                    }
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < apt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count >= 110) //(l_count > k_count)
            {
                if(l_count > k_count) k_count = l_count;
                // knum = 0;
                key_guess[knum] = k;
                knum++;
                // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fp = fopen("Result_SPECK32.txt", "a");
                // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fclose(fp);
            }
            /*
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fp = fopen("Result_SPECK32.txt", "a");
                // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                // fclose(fp);
            }
            */
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < apt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "Verified key: %llx, encoding count: %d\n", key_guess[k], k_count);
        printf("Verified key: %llx, encoding count: %d\n", key_guess[k], k_count);
        reduced_key[k] = key_guess[k];
        key_count++;
    }
    *reduced_key_count = key_count;
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "Verified key count: %d\n", key_count);
    printf("Verified key count: %d\n", key_count);
    fclose(fp); 
}

void DCA128_KeyCheck(u64 *key_candidate, int key_candidate_count, Aff128 Aff)
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, j, s, ts, l, b, bit, i;
    u64 k, kkey;
    V128 state1;
    V128 state2;
    
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key64, roundkey);

    int key_count = 0;

    int f = 0, g = 0;
    u64 map[pt][2]; 
    FILE *fp = fopen("Result_SPECK128.txt", "a");
    printf("------Final Key Verification-----\n");  
    fprintf(fp, "------Final Key Verification-----\n"); 

    u64 key_guess[200] = {0};
    int k_count = 0;
    double k_max = 0.0;
    int knum = 0;
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < pt; x++) // adaptive inputs
        {
            x0 = (u64)(x) << 5;
            x1 = 1;
            
            y1 = ROR128(ROL128(x0 - x1, 8) ^ x1, 3);
            y0 = ROL128((ROL128(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL128(y1, 3) ^ (ROR128(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR128((ROR128(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state1.V[0] = z0; 
            state1.V[1] = z1;
            MatMulVecM128(Aff.Mat, state1, &state2);
            map[x][0] = state2.V[0] ^ Aff.Vec.V[0];//add
            map[x][1] = state2.V[1] ^ Aff.Vec.V[1];
        }
        double k_score = 0.0;
        double L_score[pt] = {0.0};
        int L_count[pt] = {0};
        int l_count = 0;
        double l_max = 0.0;
        for(l = 1; l < pt; l++)
        {  
            double score = 0.0;
            double j_max = 0.0;
            int j_count = 0;
            for(j = 0; j < tb128; j++) // samples of traces
            {
                int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                for(x = 0; x < pt; x++)
                { 
                    if(j < 64)
                    {
                        if(map[x][0] & idM64[j]) ts = 1; // traces
                        else ts = 0;
                    }
                    else
                    {
                        if(map[x][1] & idM64[j - 64]) ts = 1; // traces
                        else ts = 0;
                    }
                    ////// correlation
                    f = ts;
                    if(f) Nf1++;
                    else Nf0++;
                    
                    g = xor[l & x];
                    if(g) Ng1++;
                    else Ng0++;
                    
                    if(f == 1 && g == 1) N11++;
                    else if(f == 0 & g == 0) N00++;
                    else if(f == 1 & g == 0) N10++;
                    else N01++;
                }
                if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                else score = 0.0;
                if(score > j_max)
                {
                    j_count = 0;
                    j_max = score;
                    j_count++;
                }
                else if(fabs(score - j_max) < EPS)
                {
                    j_count++;
                }
            }
            L_score[l] = j_max;
            L_count[l] = j_count;
            if(L_score[l] > l_max) l_max = L_score[l];
        }
        k_score = l_max;
        if(fabs(k_score - k_max) < EPS)
        {
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            if(l_count > k_count)
            {
                k_count = l_count;
                knum = 0;
                key_guess[knum] = k;
                knum++;
                printf("key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fprintf(fp, "key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
            else if(l_count == k_count)
            {
                key_guess[knum] = k;
                knum++;
                printf("key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                fprintf(fp, "key guess: %llx, encoding count: %d, correlation: %f\n", k, l_count, k_score);
            }
        }
        else if(k_score > k_max) 
        {
            k_max = k_score; 
            knum = 0;
            for(l = 1; l < pt; l++)
            {
                if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                {
                    l_count += L_count[l];
                }
            }
            k_count = l_count;
            key_guess[knum] = k;
            knum++;
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "simulation WB-SPECK128 recoverd key: %llx, encoding count: %d\n", key_guess[k], k_count);
        printf("simulation WB-SPECK128 recoverd key: %llx, encoding count: %d\n", key_guess[k], k_count);
        key_count++;
    }
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "simulation WB-SPECK128 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK128 recoverd key count: %d\n", key_count);
    fclose(fp); 
}

void DCA128()
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, i, j, s, ts, l, b, bit, q, len;
    u64 k, kkey;
    V128 state1;
    V128 state2;
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key64, roundkey);
    ////

    int f = 0, g = 0;
    Aff128 Aff, Aff_inv;
    genaffinepairM128(&Aff, &Aff_inv);
    u64 map[pt][2]; 

    FILE *fp; 

    u64 key_candidate[200] = {0};
    int key_candidate_count = 1;
    u64 temp_key_candidate[200] = {0};
    int temp_key_candidate_count = 0;
    for(q = 0; q < 16; q++)
    {
        memset(temp_key_candidate, 0, sizeof(temp_key_candidate));
        temp_key_candidate_count = 0;
        for(s = 0; s < key_candidate_count; s++)
        {
            u64 key_g = key_candidate[s];
            u64 key_guess[20] = {0};
            int k_count = 0;
            double k_max = 0.0;
            int knum = 0;
            for(kkey = ktb; kkey <= kte; kkey++) // key
            {
                k = (kkey << 4 * q) ^ key_g;
                for(x = 0; x < pt; x++) // adaptive inputs
                {
                    x0 = (u64)(x) << 4 * q;
                    x1 = 0;
                    
                    y1 = ROR128(ROL128(x0 - x1, 8) ^ x1, 3);
                    y0 = ROL128((ROL128(x0 - x1, 8) ^ k) - y1, 8);
                    
                    z1 = ROL128(y1, 3) ^ (ROR128(y0, 8) + y1) ^ roundkey[0];
                    z0 = ROR128((ROR128(y0, 8) + y1) ^ roundkey[0], 8) + z1;

                    state1.V[0] = z0; 
                    state1.V[1] = z1;
                    MatMulVecM128(Aff.Mat, state1, &state2);
                    map[x][0] = state2.V[0] ^ Aff.Vec.V[0];//add
                    map[x][1] = state2.V[1] ^ Aff.Vec.V[1];
                }
                double k_score = 0.0;
                double L_score[pt] = {0.0};
                int L_count[pt] = {0};
                int l_count = 0;
                double l_max = 0.0;
                for(l = 1; l < pt; l++)
                {  
                    double score = 0.0;
                    double j_max = 0.0;
                    int j_count = 0;
                    for(j = 0; j < tb128; j++) // samples of traces
                    {
                        int Nf0 = 0, Nf1 = 0, Ng0 = 0, Ng1 = 0, N00 = 0, N01 = 0, N10 = 0, N11 = 0;
                        for(x = 0; x < pt; x++)
                        { 
                            if(j < 64)
                            {
                                if(map[x][0] & idM64[j]) ts = 1; // traces
                                else ts = 0;
                            }
                            else
                            {
                                if(map[x][1] & idM64[j - 64]) ts = 1; // traces
                                else ts = 0;
                            }
                            ////// correlation
                            f = ts;
                            if(f) Nf1++;
                            else Nf0++;
                            
                            g = xor[l & x];
                            if(g) Ng1++;
                            else Ng0++;
                            
                            if(f == 1 && g == 1) N11++;
                            else if(f == 0 & g == 0) N00++;
                            else if(f == 1 & g == 0) N10++;
                            else N01++;
                        }
                        if(Nf1 && Nf0 && Ng1 && Ng0) score = abs((N11 * N00 - N10 * N01)) * 1.0 / (sqrt(Nf1) * sqrt(Nf0) * sqrt(Ng1) * sqrt(Ng0));
                        else score = 0.0;
                        if(score > j_max)
                        {
                            j_count = 0;
                            j_max = score;
                            j_count++;
                        }
                        else if(fabs(score - j_max) < EPS)
                        {
                            j_count++;
                        }
                    }
                    L_score[l] = j_max;
                    L_count[l] = j_count;
                    if(L_score[l] > l_max) l_max = L_score[l];
                }
                k_score = l_max;
                if(fabs(k_score - k_max) < EPS)
                {
                    for(l = 1; l < pt; l++)
                    {
                        if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                        {
                            l_count += L_count[l];
                        }
                    }
                    if(l_count >= 110) //(l_count > k_count)
                    {
                        if(l_count > k_count) k_count = l_count;
                        // knum = 0;
                        key_guess[knum] = k;
                        knum++;
                        // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fp = fopen("Result_SPECK32.txt", "a");
                        // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fclose(fp);
                    }
                    /*
                    else if(l_count == k_count)
                    {
                        key_guess[knum] = k;
                        knum++;
                        // printf("key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fp = fopen("Result_SPECK32.txt", "a");
                        // fprintf(fp, "key guess: %.2x, encoding count: %d, correlation: %f\n", k, l_count, k_score);
                        // fclose(fp);
                    }
                    */
                }
                else if(k_score > k_max) 
                {
                    k_max = k_score; 
                    knum = 0;
                    for(l = 1; l < pt; l++)
                    {
                        if(fabs(L_score[l] - l_max) < EPS) // each encoding: key guess
                        {
                            l_count += L_count[l];
                        }
                    }
                    k_count = l_count;
                    key_guess[knum] = k;
                    knum++;
                }
            }
            fp = fopen("Result_SPECK128.txt", "a");
            printf("------Partial Key Guess-----\n");  
            fprintf(fp, "------Partial Key Guess-----\n");
            printf("------\n");  
            fprintf(fp, "------\n");  
            for(k = 0; k < knum; k++)
            {
                fprintf(fp, "No.%d nibble: recoverd key = %llx, encoding count = %d, correlation: %f\n", q, key_guess[k], k_count, k_max);
                printf("No.%d nibble: recoverd key = %llx, encoding count = %d, correlation: %f\n", q, key_guess[k], k_count, k_max);
                temp_key_candidate[temp_key_candidate_count] = key_guess[k];
                temp_key_candidate_count++;
            }
            printf("------\n"); 
            fprintf(fp, "------\n");
            fclose(fp);  
        }
        len = sizeof(temp_key_candidate) / sizeof(temp_key_candidate[0]);
        memcpy(key_candidate, temp_key_candidate, len * sizeof(u64));
        key_candidate_count = temp_key_candidate_count;
        if(q > 0) sub_DCA128_KeyCheck(key_candidate, key_candidate_count, Aff, q, key_candidate, &key_candidate_count);
    }
    fp = fopen("Result_SPECK128.txt", "a");
    printf("------Possible Key-----\n");  
    fprintf(fp, "------Possible Key-----\n");  
    for(i = 0; i < key_candidate_count; i++)
    {
        printf("%llx\n", key_candidate[i]);
        fprintf(fp, "%llx\n", key_candidate[i]);
    }
    printf("Count: %d\n\n", key_candidate_count);  
    fprintf(fp, "Count: %d\n\n", key_candidate_count);

    fclose(fp); 

    DCA128_KeyCheck(key_candidate, key_candidate_count, Aff);
}