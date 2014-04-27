; ------------------------------------------------------------------------------------------------
; Creation date: 2014.04.16
; Creators: Youcef Lemsafer
; Authors: Youcef Lemsafer
; ------------------------------------------------------------------------------------------------

; ------------------------------------------------------------------------------
bits 64
default rel
; ------------------------------------------------------------------------------

extern printf


; ------------------------------------------------------------------------------------------------
SECTION .data
; ------------------------------------------------------------------------------------------------
printfmt:   db   '%llu', 0dh, 0ah, 0
odd_primes: dq   3,  5,  7,  11, 13, 17, 19, 23, \
                29, 31, 37, 41, 43, 47, 53, 59, \
                61, 67, 71, 73, 79, 83, 89, 97, \
                0
; L = 97, trial division limit
td_bound:  dq    97*97

%macro  findd    2                  ; %1 is d
                                    ; %2 is p
    mov     rdx, %2
    mov     rcx, %1
    call    jacobi_symbol
    test    rax, rax                ; if jacobi_symbol(d, p) = 0 then we've found a factor
    jz      .the_end
    mov     r11, %1
    test    rax, rax
    js      .found_d
%endmacro

; ------------------------------------------------------------------------------------------------
SECTION .code
; ------------------------------------------------------------------------------------------------

global is_prime
global jacobi_symbol
global is_slprp

; ------------------------------------------------------------------------------------------------
; bool is_prime( unsigned long long p )
; ------------------------------------------------------------------------------------------------
is_prime:
; ------------------------------------------------------------------------------------------------
    sub     rsp, 10h
    ; At this point p is in rcx
    ; If p is 2 return true
    mov     rdx, rcx
    xor     rdx, 2
    jz      .ret_is_prime_1
    ; If p != 1 (mod 2) or (p == 1) return false
    mov     rdx, rcx
    and     rdx, 1
    jz      .return_is_prime_0      ; p != 1 (mod 2)
    mov     rdx, rcx
    xor     rdx, 1
    jz      .return_is_prime_0      ; p == 1 (mod 2)

    ; OK, at this point p is odd and >= 3
    ; (0) Step zero: trial division by odd primes < L
    lea     r8, [odd_primes]
.td_loop:
    mov     r9, [r8]
    cmp     r9, 0
    je      .end_td_loop
    cmp     rcx, r9                 ; before dividing we test whether p is a prime <= L
    je      .ret_is_prime_1
    mov     rax, rcx                ; perform division. rax <- dividend, divisor already in r9
    xor     rdx, rdx                ; before performing division we have to reset rdx
    div     r9                      ; rdx <- remainder
    cmp     rdx, 0
    je      .return_is_prime_0      ; trial division succeeded, the number is not a prime

    add     r8, 8                   ; 8 = sizeof(odd_primes[0]), move to next item in odd_primes
    jmp     .td_loop

.end_td_loop:
    cmp     rcx, [td_bound]         ; at this point if p < L^2 we know for sure that it is a prime
    jb      .ret_is_prime_1

    ; (1) Step one: if p is not 2-sprp it is composite so the test fails
    call    is_2sprp                ; @todo: is it possible to avoid this call?
    test    al, 1
    jz      .the_end

    ; (2) Step two: search for D in {5, -7, 9, -11, ... } such that jacobi_symbol(D, p) = -1
    mov     [rsp + 8], rcx
    findd    5, rcx
    findd   -7, [rsp + 8]
    findd    9, [rsp + 8]
    findd  -11, [rsp + 8]
    findd   13, [rsp + 8]
    findd  -15, [rsp + 8]
    findd   17, [rsp + 8]
    findd  -19, [rsp + 8]
    findd   21, [rsp + 8]
    findd  -23, [rsp + 8]
    findd   25, [rsp + 8]
    findd  -27, [rsp + 8]
    findd   29, [rsp + 8]

    mov     r11, -31                ; we rely on the fact that jacobi_symbol does not use r11
.loop_find_d:
    mov     rcx, r11                ; rcx <- d
    mov     rdx, [rsp + 8]          ; rdx <- p
    call    jacobi_symbol
    test    rax, rax
    jz      .the_end                ; (D|p) = 0 => p is composite
    test    rax, rax
    js      .found_d
    test    r11, r11
    js      .d_is_negative
    add     r11, 2
    jmp     .negate_d
.d_is_negative:
    sub     r11, 2
.negate_d:
    neg     r11
    jmp     .loop_find_d

.found_d:
    ; Last step: if p is a strong Lucas probable prime for P=1 and Q=(1-D)/4
    ; then it is a prime
    mov     rcx, [rsp + 8]
    mov     rdx, r11
    call    is_slprp
    jmp     .the_end

.return_is_prime_0:
    mov     al, 0
    jmp     .the_end
.ret_is_prime_1:
    mov     al, 1
.the_end:
    add     rsp, 10h
    ret
; ------------------------------------------------------------------------------------------------
; end is_prime
; ------------------------------------------------------------------------------------------------


; ------------------------------------------------------------------------------------------------
; Returns 1 if and only if p is a base 2 strong probable prime
; otherwise returns 0.
; Important: p must be odd >= 3
; bool is_2sprp(unsigned long long p)
; return value in al
; p is expected in rcx
; ------------------------------------------------------------------------------------------------
is_2sprp:
; ------------------------------------------------------------------------------------------------
    ; Compute r and s such that p - 1 = (2^r).s where s is odd,
    ; at the end of the computation r is in r8 and s is in r9.
    mov     r9, rcx
    dec     r9                      ; r9 = s <- p - 1
    shr     r9, 1                   ; r9 = s <- s / 2   |
    mov     r8, 1                   ; r8 = r <- 1       |-> because we know for sure that p-1 is even so r >= 1
.div_by_2_loop:
    mov     r10, r9                 ; if r9 is odd then end the loop, otherwise...
    and     r10, 1
    jnz     .end_div_by_2_loop
    shr     r9, 1                   ; s <- s / 2
    inc     r8                      ; ++r
    jmp     .div_by_2_loop

.end_div_by_2_loop:
    ; now s is in r9 and r is in r8...
    ; p is 2-sprp if and only if 2^s = 1 (mod p) or 2^(s.2^i) = -1 (mod p)
    ; for some i = 0, 1, ..., r-1
    ; To compute 2^s (mod p) we use binary ladder exponentiation (left-right
    ; form) from [CP](§9.3.1).
    mov     r10, rcx                ; backup rcx value in r10
    bsr     rcx, r9                 ; put index of s's high bit in rcx
    mov     rax, 1
    shl     rax, cl
    mov     rcx, r10                ; restore rcx
    mov     r10, rax
    shr     r10, 1                  ; because the loop starts with j = D-2
    mov     rax, 2

.ble_loop:                          ; binary ladder exponentiation loop
    test    r10, r10
    jz      .end_ble_loop
    mul     rax                     ; rax <- rax * rax
    div     rcx
    mov     rax, rdx                ; rax <- rax mod p
    test    r10, r9
    jz      .current_bit_not_set
    shl     rax, 1                  ; rax <- 2 * rax
    jc      .mul2_overflow
    xor     rdx, rdx
    jmp     .no_mul2_overflow
.mul2_overflow:
    mov     rdx, 1
.no_mul2_overflow:
    div     rcx
    mov     rax, rdx                ; rax <- rax mod p
.current_bit_not_set:
    shr     r10, 1
    jmp     .ble_loop
.end_ble_loop:
    cmp     rax, 1
    je      .ret_instruction        ; since rax == 1 we can leace immediately
    mov     r10, rcx
    dec     r10                     ; r10 <- p-1
.sprp_loop:
    cmp     r10, rax
    je      .return_true
    dec     r8
    jz      .end_sprp_loop
    mul     rax                     ; rax <- rax * rax
    div     rcx
    mov     rax, rdx                ; rax <- rax mod p
    jmp     .sprp_loop
.end_sprp_loop:
    mov     al, 0
    jmp     .ret_instruction
.return_true:
    mov     al, 1
.ret_instruction:
    ret
; ------------------------------------------------------------------------------------------------
; end is_2sprp
; ------------------------------------------------------------------------------------------------


; ------------------------------------------------------------------------------------------------
; Computes the Jacobi symbol (a|m)
; for odd positive integer m and integer a
; ------------------------------------------------------------------------------------------------
; int64 jacobi_symbol(int64 a, uint64 m)
; a in rcx, m in rdx, return value in rax
; used registers: rax, rcx, rdx, r8, r9, r10
; ------------------------------------------------------------------------------------------------
; Note: To ease implementation, when a < 0 I call ujacobi_symbol
; with negated a, then I use the following identity:
;   (-x|m) = (-1)^((m-1)/2).(x|m)
; ------------------------------------------------------------------------------------------------
jacobi_symbol:
; ------------------------------------------------------------------------------------------------
    test    rcx, rcx
    js      .negative_a
    call    ujacobi_symbol
    ret
.negative_a:
    neg     rcx
    mov     r10, rdx                ; backup m in r10 (it's OK because ujacobi_symbol doesn't use
                                    ; that register)
    call    ujacobi_symbol
    dec     r10
    shr     r10, 1
    test    r10, 1
    jz      .do_not_neg_rax
    neg     rax
.do_not_neg_rax:
    ret
; ------------------------------------------------------------------------------------------------
; end jacobi_symbol
; ------------------------------------------------------------------------------------------------


; ------------------------------------------------------------------------------------------------
; Computes the Jacobi symbol (a|m)
; for odd positive integer m and integer a >= 0
; ------------------------------------------------------------------------------------------------
; int64 ujacobi_symbol(uint64 a, uint64 m)
; a in rcx, m in rdx, return value in rax
; used registers: rax, rcx, rdx, r8, r9
; ------------------------------------------------------------------------------------------------
ujacobi_symbol:
; ------------------------------------------------------------------------------------------------
    mov     r8, rdx                 ; backup m in r8
    xor     rdx, rdx
    mov     rax, rcx
    div     r8                      ; a mod m is in rdx
    mov     rcx, rdx                ; a <- a mod m
    mov     rdx, r8                 ; restore m
    mov     r8, 1                   ; t = 1
.loop_a_not_zero:
    test    rcx, rcx
    jz      .end_loop_a_not_zero
.loop_a_even:
    test    rcx, 1
    jnz     .end_loop_a_even
    shr     rcx, 1                  ; a <- a / 2
    mov     r9, rdx                 ; if m is 3 or 5 modulo 8 negate t
    and     r9, 7
    cmp     r9, 3
    je      .neg_t
    cmp     r9, 5
    je      .neg_t
    jmp     .loop_a_even
.neg_t:
    neg     r8
    jmp     .loop_a_even
.end_loop_a_even:
    xor     rcx, rdx                ;|
    xor     rdx, rcx                ;| -> swap a and m
    xor     rcx, rdx                ;|
    mov     r9, rcx                 ; if both a and m are 3 modulo 4 negate t
    and     r9, 3
    cmp     r9, 3
    jne     .not_3_mod_4
    mov     r9, rdx
    and     r9, 3
    cmp     r9, 3
    jne     .not_3_mod_4
    neg     r8
.not_3_mod_4:
    mov     r9, rdx                 ; a <- a mod m
    xor     rdx, rdx
    mov     rax, rcx
    div     r9
    mov     rcx, rdx
    mov     rdx, r9
    jmp     .loop_a_not_zero
.end_loop_a_not_zero:
    cmp     rdx, 1
    je      .return_t
    mov     rax, 0
    ret
.return_t:
    mov     rax, r8
    ret 
; ------------------------------------------------------------------------------------------------
; end ujacobi_symbol
; ------------------------------------------------------------------------------------------------


; ------------------------------------------------------------------------------
; We assume that n is odd > 1 and that Q = (1 - D)/4
; where D is the first element in {5, -7, 9, -11, ...}
; such that the Jacobi symbol (D/n) = -1.
; ------------------------------------------------------------------------------
; bool is_slprp(uint64 n, int64 D);
; ------------------------------------------------------------------------------
; Recall that:
;     - on Linux/Unix n is in rdi, D is in rsi
;     - on Windows n is in rcx and D is in rdx
; Implementation:
;   (1) Compute Q = (1 - D)/4
;   (2) Factor n + 1 into the form d.2^s
;       Use the following relations:
;           U(2k) = U(k).V(k)
;           V(2k) = V(k)^2 - 2.Q^k
;           U(2k+1) = (U(2k) + V(2k)) / 2
;           V(2k+1) = (D.U(2k) + V(2k)) / 2
;       to compute U(d) mod n.
;       n is slprp if and only if one of the following conditions holds:
;           U(d) = 0 (mod n)
;           V(d.2^r) = 0 (mod n) for some r < s
; ------------------------------------------------------------------------------
is_slprp:
    sub     rsp, 40h
%ifndef WIN64
    mov     rcx, rdi    ; this first two lines is to make it look like Windows,
    mov     rdx, rsi    ; they have to be removed if we are on Windows
%endif
    mov     [rsp], rcx
    mov     [rsp + 8h], rdx
    mov     rax, 1                      ; compute and store 2^63 used below
    shl     rax, 63                     ; for division by 2 of unsigned 128-bit
    mov     qword [rsp + 28h], rax      ; number rdx:rax (because we are not
                                        ; allowed to use immediate 1 << 63)
    mov     rax, rdx
    neg     rax
    add     rax, 1
    sar     rax, 2      ; after that Q is in rax
    ; Since we are working modulo n here we can make Q positive by adding n
    ; (we assume that |Q| < n).
    test    rax, rax
    jns     .q_is_not_negative
    add     rax, rcx
.q_is_not_negative:
    mov     [rsp + 10h], rax    ; rsp + 10h <- Q

    ; Same for D: if D is < 0 then D <- D + n.
    test    rdx, rdx
    jns     .D_is_not_negative
    add     rdx, rcx
.D_is_not_negative:
    mov     qword [rsp + 8h], rdx       ; rsp + 8h <- D


    ; Factor n + 1 into d.2^s, d in r8 and s in r9
    mov     r10, 8000000000000000h
    mov     r8, rcx     ; initialize d = n + 1
    add     r8, 1       ; we do not test for overflow here since 2^64-1 =
                        ; 3.5.17.257.641.65537.P7 cannot end here (because
                        ; we assume that trial division by small primes has
                        ; been performed).
    mov     r9, 0       ; initialize s = 0
.div2:
    test    r8, 1
    jnz     .div2_end
    shr     r8, 1
    shr     r10, 1
    inc     r9
    jmp     .div2
.div2_end:                           ; at this point d is in r8 and s is in r9

    ; Compute U(d) mod n and V(d) mod n

.adj_d_msbmask:
    test    r8, r10
    jnz     .adj_d_msbmask_end
    shr     r10, 1
    jmp     .adj_d_msbmask
.adj_d_msbmask_end:

    ; initialization
    mov     rcx, 0                      ; U = U(0) = 0
    mov     r11, 2                      ; V = V(0) = 2
    mov     qword [rsp + 18h], 1        ; q is the k-th power of Q
.loop_d_msbmask:                        ; loop until r10 is 0
    test    r10, r10
    jz      .loop_d_msbmask_end
    mov     rax, rcx                    ; compute UV
    mul     r11
    div     qword [rsp]
    mov     rcx, rdx                    ; U = UV mod n

    mov     rax, r11                    ; compute V^2 - 2Q^k = V^2 - 2q (mod n)
    mul     r11
    div     qword [rsp]
    mov     r11, rdx                    ; r11 <- V^2 (mod n)

    cmp     r11, qword [rsp + 18h]      ; if (V^2 < q) then V^2 - q would be < 0
    jb      .v2_b_q
    jae     .v2_ae_q
.v2_b_q:
    mov     rdx, qword [rsp + 18h]
    sub     rdx, r11
    mov     r11, qword [rsp]
    sub     r11, rdx
    jmp     .after_v2_ae_q
.v2_ae_q:
    sub     r11, qword [rsp + 18h]
.after_v2_ae_q:

    cmp     r11, qword [rsp + 18h]      ; if (V^2 - q < q) then V^2 - 2q would
                                        ; be negative
    jb      .v2mq_b_q
    jae     .v2mq_ae_q
.v2mq_b_q:
    mov     rdx, qword [rsp + 18h]
    sub     rdx, r11
    mov     r11, qword [rsp]
    sub     r11, rdx
    jmp     .after_v2mq_ae_q
.v2mq_ae_q:
    sub     r11, qword [rsp + 18h]
.after_v2mq_ae_q:                       ; now r11 contains V^2 - 2q (mod n)


    mov     rax, qword [rsp + 18h]      ; update q for next round:
    mul     rax
    div     qword [rsp]
    mov     qword [rsp + 18h], rdx      ; q <- q^2 (mod n)

    test    r10, r8
    jz      .d_current_bit_not_set

    ; U(2k+1) = (U(2k) + V(2k)) / 2
    mov     rax, rcx
    xor     rdx, rdx
    add     rax, r11
    adc     rdx, rdx
    div     qword [rsp]
    mov     rax, rdx
    xor     rdx, rdx
    test    rax, 1
    jz      .uv_sum_is_even
    add     rax, qword [rsp]
    adc     rdx, rdx
.uv_sum_is_even:
    shr     rax, 1                      ; division by 2 of the unsigned 128-bit
    shr     rdx, 1                      ; integer rdx:rax
    jnc     .no_carry_1
    xor     rax, qword [rsp + 28h]
.no_carry_1:
    div     qword [rsp]
    mov     qword [rsp + 20h], rdx      ; backup U(2k+1), we still need U(2k)

    ; V(2k+1) = (D.U(2k) + V(2k)) / 2
    mov     rax, rcx                    ; rax <- U(2k)
    mul     qword [rsp + 8h]            ; rax <- U(2k) * D
    div     qword [rsp]
    mov     rax, rdx
    xor     rdx, rdx
    add     rax, r11
    adc     rdx, rdx
    div     qword [rsp]
    mov     rax, rdx
    xor     rdx, rdx
    test    rax, 1
    jz      .duv_sum_is_even
    add     rax, qword [rsp]
    adc     rdx, rdx
.duv_sum_is_even:
    shr     rax, 1                      ; division by 2 of the unsigned 128-bit
    shr     rdx, 1                      ; integer rdx:rax
    jnc     .no_carry_2
    xor     rax, qword [rsp + 28h]
.no_carry_2:
    div     qword [rsp]
    mov     r11, rdx                    ; r11 <- V(2k+1)
    mov     rcx, qword [rsp + 20h]      ; rcx <- U(2k+1)

    ; q <- q * Q (mod n)
    mov     rax, qword [rsp + 18h]      ; rax <- q
    mul     qword [rsp + 10h]           ; multiply by Q
    div     qword [rsp]
    mov     qword [rsp + 18h], rdx      ; q <- Q^k (mod n)

.d_current_bit_not_set:
    shr     r10, 1
    jmp     .loop_d_msbmask

.loop_d_msbmask_end:

    ; Now U(d) is in rcx and V(d) is in r11 and Q^d is in [rsp + 18h]
    test    rcx, rcx                    ; U(d) = 0 (mod n) => n is slprp
    jz      .is_slprp_true
    test    r11, r11                    ; V = 0 (mod n) => n is slprp
    jz      .is_slprp_true

    dec     r9                          ; Loop for i in {1, ..., s-1}
.loop_on_s:
    test    r9, r9
    jz      .loop_on_s_end
    mov     rax, r11                    ; *begin* V <- V^2 - 2q (mod n)
    mul     rax                         ; rdx:rax <- V^2
    div     qword [rsp]
    mov     rax, rdx                    ; rax <- V^2 (mod n)

    ; subtract q from V^2 twice avoiding negative values
    cmp     rax, [rsp + 18h]
    jb      .v2_less_than_q
    sub     rax, [rsp + 18h]
    jmp     .after_v2_less_than_q
.v2_less_than_q:
    mov     rdx, [rsp + 18h]            ; rdx <- q
    sub     rdx, rax                    ; rdx <- rdx - V^2
    mov     rax, qword [rsp]            ; rax <- n
    sub     rax, rdx                    ; rax <- rax - rdx = n - ( q - V^2 )
.after_v2_less_than_q:
    cmp     rax, [rsp + 18h]            ; compare V^2 - q and q
    jb      .v2mq_less_than_q
    sub     rax, [rsp + 18h]            ; rax <- rax - q
    jmp     .after_v2mq_less_than_q
.v2mq_less_than_q:
    mov     rdx, [rsp + 18h]            ; rdx <- q
    sub     rdx, rax                    ; rdx <- rdx - V^2
    mov     rax, qword [rsp]            ; rax <- n
    sub     rax, rdx                    ; rax <- rax - rdx = n - ( q - V^2 )
.after_v2mq_less_than_q:
    
    test    rax, rax                    ; test rdx before moving it to r11
    jz      .is_slprp_true
    mov     r11, rax                    ; *end* V <- V^2 - 2q (mod n)

    mov     rax, [rsp + 18h]            ; *begin* q <- q^2 (mod n)
    mul     rax
    div     qword [rsp]
    mov     [rsp + 18h], rdx            ; *end* q <- q^2 (mod n)

    dec     r9
    jmp     .loop_on_s
.loop_on_s_end:
    xor     rax, rax
    jmp     .is_slprp_end

.is_slprp_true:
    mov     rax, 1

.is_slprp_end:
    add     rsp, 40h
    ret


; [CP] Richard Crandall, Carl Pomerance - Prime numbers, a computational perspective.
;       Second edition Springer 2005
; Algorithm 9.3.1 (binary ladder exponentiation (left-right form))
; Computes x^y. Assume (y_0,..., y_D-1) is the binary expansion of y > 0 where
; y_D-1 = 1 is the high bit.
;   // Initialize
;       z = x
;   // Loop over bits of y, starting with next-to-highest
;       for(D-2 >= j >= 0) {
;           z = z^2                 // For modular arithmetic, do mod N here
;           if(y_j == 1) z = zx     // For modular arithmetic, do mod N here
;       }
;       return z
