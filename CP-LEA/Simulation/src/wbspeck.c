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

void sub_LEA32_KeyCheck(u16 *key_candidate, int key_candidate_count, Aff32 Aff, int q, u16 *reduced_key, int *reduced_key_count)
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u32 state;
    
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key16, roundkey);

    int key_count = 0;

    u32 map[qpt]; 
    FILE *fp = fopen("Result_SPECK32.txt", "a");
    printf("------Partial Key Verification-----\n");  
    fprintf(fp, "------Partial Key Verification-----\n"); 

    u8 Input[qpt] = {0x0, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, 0xff};
    uint8_t Htrace[qpt][9];
    uint8_t temp;
    uint8_t trail[50][3];// Gaussian trail
    for(x = 0; x < qpt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 9; bit++)
        {
            if(Input[x] & idM8[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 40
    //  Gaussian Elimination
    for(i = 0; i < 9; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < qpt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 9; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < qpt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 9; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < qpt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 9; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = qpt;
    for(i = qpt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 9; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u16 key_guess[100] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[qpt];
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < qpt; x++) // adaptive inputs
        {
            x0 = (u16)(Input[x]) << ((q - 1) * 4);
            x1 = 0;
            
            y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
            y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
            
            z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
            z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

            state = (z0 << 16) | z1;
            map[x] = state;
            map[x] = affineU32(Aff, state); // affine encoding
        }
        int k_count = 0;
        for(j = 0; j < tb32; j++) // samples of traces
        {
            for(x = 0; x < qpt; x++)
            { 
                if(map[x] & idM32[j]) vector[x] = 1; // traces
                else vector[x] = 0;
            }
            //  Gaussian Elimination
            for(i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = qpt;
            for(i = qpt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            // fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            // printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            // fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            // printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "Verified key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        printf("Verified key: %.2x, encoding count: %d\n", key_guess[k], k_max);
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

void LEA32_KeyCheck(u16 *key_candidate, int key_candidate_count, Aff32 Aff)
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u32 state;
    
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key16, roundkey);

    int key_count = 0;

    u32 map[ppt]; 
    FILE *fp = fopen("Result_SPECK32.txt", "a");
    printf("------Final Key Verification-----\n");  
    fprintf(fp, "------Final Key Verification-----\n"); 

    u8 Input[ppt] = {0x0, 0x8, 0x4, 0x2, 0x1, 0xf};
    uint8_t Htrace[ppt][5];
    uint8_t temp;
    uint8_t trail[20][3];// Gaussian trail
    for(x = 0; x < ppt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 5; bit++)
        {
            if(Input[x] & idM4[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 4
    //  Gaussian Elimination
    for(i = 0; i < 5; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 5; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 5; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < ppt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 5; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = ppt;
    for(i = ppt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 5; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u16 key_guess[100] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[ppt];
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < ppt; x++) // adaptive inputs
        {
            x0 = (u16)(Input[x]) << 5;
            x1 = 1;
            
            y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
            y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
            
            z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
            z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

            state = (z0 << 16) | z1;
            map[x] = state;
            map[x] = affineU32(Aff, state); // affine encoding
        }
        int k_count = 0;
        for(j = 0; j < tb32; j++) // samples of traces
        {
            for(x = 0; x < ppt; x++)
            { 
                if(map[x] & idM32[j]) vector[x] = 1; // traces
                else vector[x] = 0;
            }
            //  Gaussian Elimination
            for(i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = ppt;
            for(i = ppt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "simulation WB-SPECK32 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        printf("simulation WB-SPECK32 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        key_count++;
    }
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK32 recoverd key count: %d\n", key_count);
    fclose(fp); 
}
void LEA32()
{
    u16 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, i, j, s, ts, l, b, bit, q, len;
    u32 state;
    u16 roundkey[SPECK_ROUNDS];
    Key_expand32(key16, roundkey);
    ////

    Aff32 Aff, Aff_inv;
    genaffinepairM32(&Aff, &Aff_inv);
    u32 map[ppt]; 
    
    FILE *fp; 

    u8 Input[ppt] = {0x0, 0x8, 0x4, 0x2, 0x1, 0xf};
    uint8_t Htrace[ppt][5];
    uint8_t temp;
    uint8_t trail[20][3];// Gaussian trail
    for(x = 0; x < ppt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 5; bit++)
        {
            if(Input[x] & idM4[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 4
    //  Gaussian Elimination
    for(i = 0; i < 5; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 5; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 5; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < ppt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 5; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = ppt;
    for(i = ppt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 5; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

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
            int k_max = 0;
            int knum = 0;
            uint8_t vector[ppt];
            for(kkey = ktb; kkey <= kte; kkey++) // key
            {
                k = (kkey << 4 * q) ^ key_g;
                for(x = 0; x < ppt; x++) // adaptive inputs
                {
                    x0 = (u16)(Input[x]) << 4 * q;
                    x1 = 0;
                    
                    y1 = ROR32(ROL32(x0 - x1, 7) ^ x1, 2);
                    y0 = ROL32((ROL32(x0 - x1, 7) ^ k) - y1, 7);
                    
                    z1 = ROL32(y1, 2) ^ (ROR32(y0, 7) + y1) ^ roundkey[0];
                    z0 = ROR32((ROR32(y0, 7) + y1) ^ roundkey[0], 7) + z1;

                    state = (z0 << 16) | z1;
                    map[x] = state;
                    map[x] = affineU32(Aff, state); // affine encoding
                }
                int k_count = 0;
                for(j = 0; j < tb32; j++) // samples of traces
                {
                    for(x = 0; x < ppt; x++)
                    { 
                        if(map[x] & idM32[j]) vector[x] = 1; // traces
                        else vector[x] = 0;
                    }
                    //  Gaussian Elimination
                    for(i = 0; i < Gauss_time; i++)
                    {
                        if(trail[i][0]) // addition
                        {
                            vector[trail[i][1]] ^= vector[trail[i][2]];
                        }
                        else // swap
                        {
                            temp = vector[trail[i][2]];
                            vector[trail[i][2]] = vector[trail[i][1]];
                            vector[trail[i][1]] = temp;
                        }
                    }
                    
                    // Gauss Over
                    int rAb = ppt;
                    for(i = ppt - 1; i >= 0; i--)
                    {
                        int allzero = 1;
                        if(vector[i]) allzero = 0;
                        if(allzero) rAb--;
                        else break;
                    }
                    if(rA >= rAb) // has a solusion
                    {
                        k_count++;
                    }
                }
                if(k_count && (k_count == k_max))
                {
                    key_guess[knum] = k; // kkey
                    knum++;
                    // fprintf(fp, "key guess: %x, encoding count: %d\n", kkey, k_count);
                    // printf("key guess: %x, encoding count: %d\n", kkey, k_count);
                }
                else if(k_count && (k_count > k_max))
                {
                    k_max = k_count; 
                    knum = 0;
                    key_guess[knum] = k; // kkey
                    knum++;
                    // fprintf(fp, "key guess: %x, encoding count: %d\n", kkey, k_count);
                    // printf("key guess: %x, encoding count: %d\n", kkey, k_count);
                }
            }
            fp = fopen("Result_SPECK32.txt", "a");
            printf("------Partial Key Guess-----\n");  
            fprintf(fp, "------Partial Key Guess-----\n");
            printf("------\n");  
            fprintf(fp, "------\n");  
            for(k = 0; k < knum; k++)
            {
                fprintf(fp, "No.%d nibble: recoverd key = %x, encoding count = %d\n", q, key_guess[k], k_max);
                printf("No.%d nibble: recoverd key = %x, encoding count = %d\n", q, key_guess[k], k_max);
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
        if(q > 0) sub_LEA32_KeyCheck(key_candidate, key_candidate_count, Aff, q, key_candidate, &key_candidate_count);    
    }
    fp = fopen("Result_SPECK32.txt", "a");
    printf("------Possible Key-----\n");  
    fprintf(fp, "------Possible Key-----\n");  
    for(i = 0; i < key_candidate_count; i++)
    {
        printf("%x\n", key_candidate[i]);
        fprintf(fp, "%x\n", key_candidate[i]);
    }
    printf("Count: %d\n\n", key_candidate_count);  
    fprintf(fp, "Count: %d\n\n", key_candidate_count); 
    fclose(fp); 

    LEA32_KeyCheck(key_candidate, key_candidate_count, Aff);
}

void sub_LEA64_KeyCheck(u32 *key_candidate, int key_candidate_count, Aff64 Aff, int q, u32 *reduced_key, int *reduced_key_count)
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u64 state;
    
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key32, roundkey);

    int key_count = 0;

    u64 map[qpt]; 
    FILE *fp = fopen("Result_SPECK64.txt", "a");
    printf("------Partial Key Verification-----\n");  
    fprintf(fp, "------Partial Key Verification-----\n"); 

    u8 Input[qpt] = {0x0, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, 0xff};
    uint8_t Htrace[qpt][9];
    uint8_t temp;
    uint8_t trail[50][3];// Gaussian trail
    for(x = 0; x < qpt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 9; bit++)
        {
            if(Input[x] & idM8[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 40
    //  Gaussian Elimination
    for(i = 0; i < 9; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < qpt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 9; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < qpt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 9; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < qpt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 9; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = qpt;
    for(i = qpt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 9; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u32 key_guess[100] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[qpt];
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < qpt; x++) // adaptive inputs
        {
            x0 = (u32)(Input[x]) << ((q - 1) * 4);
            x1 = 0;
            
            y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
            y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state = ((u64)(z0) << 32) | z1;
            map[x] = state;
            map[x] = affineU64(Aff, state); // affine encoding
        }
        int k_count = 0;
        for(j = 0; j < tb64; j++) // samples of traces
        {
            for(x = 0; x < qpt; x++)
            { 
                if(map[x] & idM64[j]) vector[x] = 1; // traces
                else vector[x] = 0;
            }
            //  Gaussian Elimination
            for(i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = qpt;
            for(i = qpt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            // fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            // printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            // fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            // printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "Verified key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        printf("Verified key: %.2x, encoding count: %d\n", key_guess[k], k_max);
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
void LEA64_KeyCheck(u32 *key_candidate, int key_candidate_count, Aff64 Aff)
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, j, s, ts, l, b, bit, i;
    u64 state;
    
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key32, roundkey);

    int key_count = 0;

    u64 map[ppt]; 
    FILE *fp = fopen("Result_SPECK64.txt", "a");
    printf("------Final Key Verification-----\n");  
    fprintf(fp, "------Final Key Verification-----\n"); 

    u8 Input[ppt] = {0x0, 0x8, 0x4, 0x2, 0x1, 0xf};
    uint8_t Htrace[ppt][5];
    uint8_t temp;
    uint8_t trail[20][3];// Gaussian trail
    for(x = 0; x < ppt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 5; bit++)
        {
            if(Input[x] & idM4[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 4
    //  Gaussian Elimination
    for(i = 0; i < 5; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 5; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 5; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < ppt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 5; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = ppt;
    for(i = ppt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 5; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u32 key_guess[100] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[ppt];
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < ppt; x++) // adaptive inputs
        {
            x0 = (u32)(Input[x]) << 5;
            x1 = 1;
            
            y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
            y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
            
            z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
            z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

            state = ((u64)(z0) << 32) | z1;
            map[x] = state;
            map[x] = affineU64(Aff, state); // affine encoding
        }
        int k_count = 0;
        for(j = 0; j < tb64; j++) // samples of traces
        {
            for(x = 0; x < ppt; x++)
            { 
                if(map[x] & idM64[j]) vector[x] = 1; // traces
                else vector[x] = 0;
            }
            //  Gaussian Elimination
            for(i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = ppt;
            for(i = ppt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            fprintf(fp, "key guess: %.2x, encoding count: %d\n", k, k_count);
            printf("key guess: %.2x, encoding count: %d\n", k, k_count);
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "simulation WB-SPECK64 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        printf("simulation WB-SPECK64 recoverd key: %.2x, encoding count: %d\n", key_guess[k], k_max);
        key_count++;
    }
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "simulation WB-SPECK64 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK64 recoverd key count: %d\n", key_count);
    fclose(fp); 
}
void LEA64()
{
    u32 x0, x1, y0, y1, z0, z1;
    int k, kkey, x, i, j, s, ts, l, b, bit, q, len;
    u64 state;
    u32 roundkey[SPECK_ROUNDS];
    Key_expand64(key32, roundkey);
    ////

    Aff64 Aff, Aff_inv;
    genaffinepairM64(&Aff, &Aff_inv);
    u64 map[ppt]; 
    
    FILE *fp;

    u8 Input[ppt] = {0x0, 0x8, 0x4, 0x2, 0x1, 0xf};
    uint8_t Htrace[ppt][5];
    uint8_t temp;
    uint8_t trail[20][3];// Gaussian trail
    for(x = 0; x < ppt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 5; bit++)
        {
            if(Input[x] & idM4[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 4
    //  Gaussian Elimination
    for(i = 0; i < 5; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 5; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 5; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < ppt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 5; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = ppt;
    for(i = ppt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 5; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

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
            int k_max = 0;
            int knum = 0;
            uint8_t vector[ppt];
            for(kkey = ktb; kkey <= kte; kkey++) // key
            {
                k = (kkey << 4 * q) ^ key_g;
                for(x = 0; x < ppt; x++) // adaptive inputs
                {
                    x0 = (u32)(Input[x]) << 4 * q;
                    x1 = 0;
                    
                    y1 = ROR64(ROL64(x0 - x1, 8) ^ x1, 3);
                    y0 = ROL64((ROL64(x0 - x1, 8) ^ k) - y1, 8);
                    
                    z1 = ROL64(y1, 3) ^ (ROR64(y0, 8) + y1) ^ roundkey[0];
                    z0 = ROR64((ROR64(y0, 8) + y1) ^ roundkey[0], 8) + z1;

                    state = ((u64)(z0) << 32) | z1;
                    map[x] = state;
                    map[x] = affineU64(Aff, state); // affine encoding
                }
                int k_count = 0;
                for(j = 0; j < tb64; j++) // samples of traces
                {
                    for(x = 0; x < ppt; x++)
                    { 
                        if(map[x] & idM64[j]) vector[x] = 1; // traces
                        else vector[x] = 0;
                    }
                    //  Gaussian Elimination
                    for(i = 0; i < Gauss_time; i++)
                    {
                        if(trail[i][0]) // addition
                        {
                            vector[trail[i][1]] ^= vector[trail[i][2]];
                        }
                        else // swap
                        {
                            temp = vector[trail[i][2]];
                            vector[trail[i][2]] = vector[trail[i][1]];
                            vector[trail[i][1]] = temp;
                        }
                    }
                    
                    // Gauss Over
                    int rAb = ppt;
                    for(i = ppt - 1; i >= 0; i--)
                    {
                        int allzero = 1;
                        if(vector[i]) allzero = 0;
                        if(allzero) rAb--;
                        else break;
                    }
                    if(rA >= rAb) // has a solusion
                    {
                        k_count++;
                    }
                }
                if(k_count && (k_count == k_max))
                {
                    key_guess[knum] = k; // kkey
                    knum++;
                    // fprintf(fp, "key guess: %x, encoding count: %d\n", kkey, k_count);
                    // printf("key guess: %x, encoding count: %d\n", kkey, k_count);
                }
                else if(k_count && (k_count > k_max))
                {
                    k_max = k_count; 
                    knum = 0;
                    key_guess[knum] = k; // kkey
                    knum++;
                    // fprintf(fp, "key guess: %x, encoding count: %d\n", kkey, k_count);
                    // printf("key guess: %x, encoding count: %d\n", kkey, k_count);
                }
            }
            fp = fopen("Result_SPECK64.txt", "a");
            printf("------Partial Key Guess-----\n");  
            fprintf(fp, "------Partial Key Guess-----\n");
            printf("------\n");  
            fprintf(fp, "------\n");  
            for(k = 0; k < knum; k++)
            {
                fprintf(fp, "No.%d nibble: recoverd key = %x, encoding count = %d\n", q, key_guess[k], k_max);
                printf("No.%d nibble: recoverd key = %x, encoding count = %d\n", q, key_guess[k], k_max);
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
        if(q > 0) sub_LEA64_KeyCheck(key_candidate, key_candidate_count, Aff, q, key_candidate, &key_candidate_count);
    }
    fp = fopen("Result_SPECK64.txt", "a");
    printf("------Possible Key-----\n");  
    fprintf(fp, "------Possible Key-----\n");  
    for(i = 0; i < key_candidate_count; i++)
    {
        printf("%x\n", key_candidate[i]);
        fprintf(fp, "%x\n", key_candidate[i]);
    }
    printf("Count: %d\n\n", key_candidate_count);  
    fprintf(fp, "Count: %d\n\n", key_candidate_count); 
    fclose(fp); 

    LEA64_KeyCheck(key_candidate, key_candidate_count, Aff);
}

void sub_LEA128_KeyCheck(u64 *key_candidate, int key_candidate_count, Aff128 Aff, int q, u64 *reduced_key, int *reduced_key_count)
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, j, s, ts, l, b, bit, i;
    u64 k, kkey;
    V128 state1;
    V128 state2;
    
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key64, roundkey);

    int key_count = 0;

    u64 map[qpt][2]; 
    FILE *fp = fopen("Result_SPECK128.txt", "a");
    printf("------Partial Key Verification-----\n");  
    fprintf(fp, "------Partial Key Verification-----\n"); 
    
    u8 Input[qpt] = {0x0, 0x80, 0x40, 0x20, 0x10, 0x8, 0x4, 0x2, 0x1, 0xff};
    uint8_t Htrace[qpt][9];
    uint8_t temp;
    uint8_t trail[50][3];// Gaussian trail
    for(x = 0; x < qpt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 9; bit++)
        {
            if(Input[x] & idM8[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 45
    //  Gaussian Elimination
    for(i = 0; i < 9; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < qpt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 9; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < qpt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 9; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < qpt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 9; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = qpt;
    for(i = qpt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 9; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u64 key_guess[200] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[qpt];
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < qpt; x++) // adaptive inputs
        {
            x0 = (u64)(Input[x]) << ((q - 1) * 4);
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
        int k_count = 0;
        for(j = 0; j < tb128; j++) // samples of traces
        {
            for(x = 0; x < qpt; x++)
            { 
                if(j < 64)
                {    
                    if(map[x][0] & idM64[j]) vector[x] = 1; // traces
                    else vector[x] = 0;
                }
                else
                {
                    if(map[x][1] & idM64[j - 64]) vector[x] = 1; // traces
                    else vector[x] = 0;
                }
            }
            //  Gaussian Elimination
            for(i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = qpt;
            for(i = qpt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            // fprintf(fp, "key guess: %llx, encoding count: %d\n", k, k_count);
            // printf("key guess: %llx, encoding count: %d\n", k, k_count);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            // fprintf(fp, "key guess: %llx, encoding count: %d\n", k, k_count);
            // printf("key guess: %llx, encoding count: %d\n", k, k_count);
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "Verified key: %llx, encoding count: %d\n", key_guess[k], k_max);
        printf("Verified key: %llx, encoding count: %d\n", key_guess[k], k_max);
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

void LEA128_KeyCheck(u64 *key_candidate, int key_candidate_count, Aff128 Aff)
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, j, s, ts, l, b, bit, i;
    u64 k, kkey;
    V128 state1;
    V128 state2;
    
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key64, roundkey);

    int key_count = 0;

    u64 map[ppt][2]; 
    FILE *fp = fopen("Result_SPECK128.txt", "a");
    printf("------Final Key Verification-----\n");  
    fprintf(fp, "------Final Key Verification-----\n"); 
    
    u8 Input[ppt] = {0x0, 0x8, 0x4, 0x2, 0x1, 0xf};
    uint8_t Htrace[ppt][5];
    uint8_t temp;
    uint8_t trail[20][3];// Gaussian trail
    for(x = 0; x < ppt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 5; bit++)
        {
            if(Input[x] & idM4[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 4
    //  Gaussian Elimination
    for(i = 0; i < 5; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 5; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 5; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < ppt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 5; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = ppt;
    for(i = ppt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 5; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

    u64 key_guess[200] = {0};
    int k_max = 0;
    int knum = 0;
    uint8_t vector[ppt];
    for(s = 0; s < key_candidate_count; s++) // key
    {
        k = key_candidate[s];
        for(x = 0; x < ppt; x++) // adaptive inputs
        {
            x0 = (u64)(Input[x]) << 5;
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
        int k_count = 0;
        for(j = 0; j < tb128; j++) // samples of traces
        {
            for(x = 0; x < ppt; x++)
            { 
                if(j < 64)
                {    
                    if(map[x][0] & idM64[j]) vector[x] = 1; // traces
                    else vector[x] = 0;
                }
                else
                {
                    if(map[x][1] & idM64[j - 64]) vector[x] = 1; // traces
                    else vector[x] = 0;
                }
            }
            //  Gaussian Elimination
            for(i = 0; i < Gauss_time; i++)
            {
                if(trail[i][0]) // addition
                {
                    vector[trail[i][1]] ^= vector[trail[i][2]];
                }
                else // swap
                {
                    temp = vector[trail[i][2]];
                    vector[trail[i][2]] = vector[trail[i][1]];
                    vector[trail[i][1]] = temp;
                }
            }
            
            // Gauss Over
            int rAb = ppt;
            for(i = ppt - 1; i >= 0; i--)
            {
                int allzero = 1;
                if(vector[i]) allzero = 0;
                if(allzero) rAb--;
                else break;
            }
            if(rA >= rAb) // has a solusion
            {
                k_count++;
            }
        }
        if(k_count && (k_count == k_max))
        {
            key_guess[knum] = k;
            knum++;
            fprintf(fp, "key guess: %llx, encoding count: %d\n", k, k_count);
            printf("key guess: %llx, encoding count: %d\n", k, k_count);
        }
        else if(k_count && (k_count > k_max))
        {
            k_max = k_count; 
            knum = 0;
            key_guess[knum] = k;
            knum++;
            fprintf(fp, "key guess: %llx, encoding count: %d\n", k, k_count);
            printf("key guess: %llx, encoding count: %d\n", k, k_count);
        }
    }
    fprintf(fp, "------\n");
    printf("------\n");
    for(k = 0; k < knum; k++)
    {
        fprintf(fp, "simulation WB-SPECK128 recoverd key: %llx, encoding count: %d\n", key_guess[k], k_max);
        printf("simulation WB-SPECK128 recoverd key: %llx, encoding count: %d\n", key_guess[k], k_max);
        key_count++;
    }
    fprintf(fp, "------\n");
    printf("------\n");

    fprintf(fp, "simulation WB-SPECK128 recoverd key count: %d\n", key_count);
    printf("simulation WB-SPECK128 recoverd key count: %d\n", key_count);
    fclose(fp); 
}

void LEA128()
{
    u64 x0, x1, y0, y1, z0, z1;
    int x, i, j, s, ts, l, b, bit, q, len;
    u64 k, kkey;
    V128 state1;
    V128 state2;
    u64 roundkey[SPECK_ROUNDS];
    Key_expand128(key64, roundkey);
    ////

    Aff128 Aff, Aff_inv;
    genaffinepairM128(&Aff, &Aff_inv);
    u64 map[ppt][2]; 

    FILE *fp; 

    u8 Input[ppt] = {0x0, 0x8, 0x4, 0x2, 0x1, 0xf};
    uint8_t Htrace[ppt][5];
    uint8_t temp;
    uint8_t trail[20][3];// Gaussian trail
    for(x = 0; x < ppt; x++)
    {
        Htrace[x][0] = 1;
        for(bit = 1; bit < 5; bit++)
        {
            if(Input[x] & idM4[bit - 1]) Htrace[x][bit] = 1;
            else Htrace[x][bit] = 0;
        }
    }
    int Gauss_time = 0; // 4
    //  Gaussian Elimination
    for(i = 0; i < 5; i++)
    {
        if(Htrace[i][i])
        {
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                { 
                    for(int r = 0; r < 5; r++)
                    {
                        Htrace[j][r] ^= Htrace[i][r];
                    }
                    trail[Gauss_time][0] = 1; //addition
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;
                }
            }
        }
        else
        {
            int flag = 0;
            for(j = i + 1; j < ppt; j++)
            {
                if(Htrace[j][i])
                {
                    for(int r = 0; r < 5; r++)
                    {
                        temp = Htrace[i][r];
                        Htrace[i][r] = Htrace[j][r];
                        Htrace[j][r] = temp;
                    }
                    trail[Gauss_time][0] = 0; //swap
                    trail[Gauss_time][1] = j;
                    trail[Gauss_time][2] = i;
                    Gauss_time++;

                    flag = 1;
                    break;
                }
            }
            if(flag)
            {
                for(j = i + 1; j < ppt; j++)
                {
                    if(Htrace[j][i])
                    { 
                        for(int r = 0; r < 5; r++)
                        {
                            Htrace[j][r] ^= Htrace[i][r];
                        }
                        trail[Gauss_time][0] = 1; //addition
                        trail[Gauss_time][1] = j;
                        trail[Gauss_time][2] = i;
                        Gauss_time++;
                    }
                }
            }
        }
    }
    // printf("Gauss_time: %d\n", Gauss_time);
    // return ;

    // the rank of A
    int rA = ppt;
    for(i = ppt - 1; i >= 0; i--)
    {
        int allzero = 1;
        for(int j = 0; j < 5; j++)
        {
            if(Htrace[i][j]) allzero = 0;
        }
        if(allzero) rA--;
        else break;
    }
    // printf("rA: %d\n", rA);
    // return ;

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
            int k_max = 0;
            int knum = 0;
            uint8_t vector[ppt];
            for(kkey = ktb; kkey <= kte; kkey++) // key
            {
                k = (kkey << 4 * q) ^ key_g;
                for(x = 0; x < ppt; x++) // adaptive inputs
                {
                    x0 = (u64)(Input[x]) << 4 * q;
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
                int k_count = 0;
                for(j = 0; j < tb128; j++) // samples of traces
                {
                    for(x = 0; x < ppt; x++)
                    { 
                        if(j < 64)
                        {    
                            if(map[x][0] & idM64[j]) vector[x] = 1; // traces
                            else vector[x] = 0;
                        }
                        else
                        {
                            if(map[x][1] & idM64[j - 64]) vector[x] = 1; // traces
                            else vector[x] = 0;
                        }
                    }
                    //  Gaussian Elimination
                    for(i = 0; i < Gauss_time; i++)
                    {
                        if(trail[i][0]) // addition
                        {
                            vector[trail[i][1]] ^= vector[trail[i][2]];
                        }
                        else // swap
                        {
                            temp = vector[trail[i][2]];
                            vector[trail[i][2]] = vector[trail[i][1]];
                            vector[trail[i][1]] = temp;
                        }
                    }
                    
                    // Gauss Over
                    int rAb = ppt;
                    for(i = ppt - 1; i >= 0; i--)
                    {
                        int allzero = 1;
                        if(vector[i]) allzero = 0;
                        if(allzero) rAb--;
                        else break;
                    }
                    if(rA >= rAb) // has a solusion
                    {
                        k_count++;
                    }
                }
                if(k_count && (k_count == k_max))
                {
                    key_guess[knum] = k; // kkey
                    knum++;
                    // fprintf(fp, "key guess: %x, encoding count: %d\n", kkey, k_count);
                    // printf("key guess: %x, encoding count: %d\n", kkey, k_count);
                }
                else if(k_count && (k_count > k_max))
                {
                    k_max = k_count; 
                    knum = 0;
                    key_guess[knum] = k; // kkey
                    knum++;
                    // fprintf(fp, "key guess: %x, encoding count: %d\n", kkey, k_count);
                    // printf("key guess: %x, encoding count: %d\n", kkey, k_count);
                }
            }
            fp = fopen("Result_SPECK128.txt", "a");
            printf("------Partial Key Guess-----\n");  
            fprintf(fp, "------Partial Key Guess-----\n");
            printf("------\n");  
            fprintf(fp, "------\n");  
            for(k = 0; k < knum; k++)
            {
                fprintf(fp, "No.%d nibble: recoverd key = %llx, encoding count = %d\n", q, key_guess[k], k_max);
                printf("No.%d nibble: recoverd key = %llx, encoding count = %d\n", q, key_guess[k], k_max);
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
        if(q > 0) sub_LEA128_KeyCheck(key_candidate, key_candidate_count, Aff, q, key_candidate, &key_candidate_count);
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

    LEA128_KeyCheck(key_candidate, key_candidate_count, Aff);
}