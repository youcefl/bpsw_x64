; ------------------------------------------------------------------------------------------------
; Creation date: 2014.04.16
; Creators: Youcef Lemsafer
; Authors: Youcef Lemsafer
; ------------------------------------------------------------------------------------------------

printf PROTO C :VARARG
is_slprp PROTO C :QWORD, :QWORD

; ------------------------------------------------------------------------------------------------
.data
; ------------------------------------------------------------------------------------------------
printfmt   db   '%llu', 0dh, 0ah, 0
odd_primes dq   3,  5,  7,  11, 13, 17, 19, 23, \
                29, 31, 37, 41, 43, 47, 53, 59, \
                61, 67, 71, 73, 79, 83, 89, 97, \
                0
; L = 97, trial division limit
td_bound  dq    97*97

fndd  macro  d, p
    mov     rdx, p
    mov     rcx, d
    call    jacobi_symbol
    test    rax, rax                ; if jacobi_symbol(D, p) = 0 then we've found a factor
    jz      the_end
    mov     r11, d
    test    rax, rax
    js      found_d
endm

; ------------------------------------------------------------------------------------------------
.code
; ------------------------------------------------------------------------------------------------

; ------------------------------------------------------------------------------------------------
; bool is_prime( unsigned long long p )
; ------------------------------------------------------------------------------------------------
is_prime proc
; ------------------------------------------------------------------------------------------------
    ; At this point p is in rcx
    ; If p is 2 return true
    mov     rdx, rcx
    xor     rdx, 2
    jz      ret_is_prime_1
    ; If p != 1 (mod 2) or (p == 1) return false
    mov     rdx, rcx
    and     rdx, 1
    jz      return_is_prime_0       ; p != 1 (mod 2)
    mov     rdx, rcx
    xor     rdx, 1
    jz      return_is_prime_0       ; p == 1 (mod 2)

    ; OK, at this point p is odd and >= 3
    ; (0) Step zero: trial division by odd primes < L
    lea     r8, odd_primes
td_loop:
    mov     r9, [r8]
    cmp     r9, 0
    je      end_td_loop
    cmp     rcx, r9                 ; before dividing we test whether p is a prime <= L
    je      ret_is_prime_1
    mov     rax, rcx                ; perform division. rax <- dividend, divisor already in r9
    xor     rdx, rdx                ; before performing division we have to reset rdx
    div     r9                      ; rdx <- remainder
    cmp     rdx, 0
    je      return_is_prime_0       ; trial division succeeded, the number is not a prime

    add     r8, 8                   ; 8 = sizeof(odd_primes[0]), move to next item in odd_primes
    jmp     td_loop

end_td_loop:
    cmp     rcx, td_bound           ; at this point if p < L^2 we know for sure that it is a prime
    jb      ret_is_prime_1

    ; (1) Step one: if p is not 2-sprp it is composite so the test fails
    call    is_2sprp                ; @todo: is it possible to avoid this call?
    test    al, 1
    jz      the_end

    ; (2) Step two: search for D in {5, -7, 9, -11, ... } such that jacobi_symbol(D, p) = -1
    mov     [rsp + 8], rcx
    fndd      5, rcx 
    fndd     -7, [rsp + 8]
    fndd      9, [rsp + 8]
    fndd    -11, [rsp + 8]
    fndd     13, [rsp + 8]
    fndd    -15, [rsp + 8]
    fndd     17, [rsp + 8]
    fndd    -19, [rsp + 8]
    fndd     21, [rsp + 8]
    fndd    -23, [rsp + 8]
    fndd     25, [rsp + 8]
    fndd    -27, [rsp + 8]
    fndd     29, [rsp + 8]
    mov     r11, -31                ; we rely on the fact that jacobi_symbol does not use r11
loop_find_d:
    mov     rcx, r11                ; rcx <- d
    mov     rdx, [rsp + 8]          ; rdx <- p
    call    jacobi_symbol
    test    rax, rax
    jz      the_end                 ; (D|p) = 0 => p is composite
    test    rax, rax
    js      found_d
    test    r11, r11
    js      d_is_negative
    add     r11, 2
    jmp     negate_d
d_is_negative:
    sub     r11, 2
negate_d:
    neg     r11
    jmp     loop_find_d

found_d:    
    mov     rcx, [rsp + 8]          ; rcx <- p
    mov     rdx, r11                ; rdx <- d
    call    is_slprp
    jmp     the_end

return_is_prime_0:
    mov     al, 0
    jmp     the_end
ret_is_prime_1:
    mov     al, 1
the_end:
    ret
; ------------------------------------------------------------------------------------------------
is_prime endp
; ------------------------------------------------------------------------------------------------


; ------------------------------------------------------------------------------------------------
; Returns 1 if and only if p is a base 2 strong probable prime
; otherwise returns 0.
; Important: p must be odd >= 3
; bool is_2sprp(unsigned long long p)
; return value in al
; p is expected in rcx
; ------------------------------------------------------------------------------------------------
is_2sprp proc
; ------------------------------------------------------------------------------------------------
    ; Compute r and s such that p - 1 = (2^r).s where s is odd,
    ; at the end of the computation r is in r8 and s is in r9.
    mov     r9, rcx
    dec     r9                      ; r9 = s <- p - 1
    shr     r9, 1                   ; r9 = s <- s / 2   |
    mov     r8, 1                   ; r8 = r <- 1       |-> because we know for sure that p-1 is even so r >= 1
div_by_2_loop:
    mov     r10, r9                 ; if r9 is odd then end the loop, otherwise...
    and     r10, 1
    jnz     end_div_by_2_loop
    shr     r9, 1                   ; s <- s / 2
    inc     r8                      ; ++r
    jmp     div_by_2_loop

end_div_by_2_loop:
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

ble_loop:                           ; binary ladder exponentiation loop
    test    r10, r10
    jz      end_ble_loop
    mul     rax                     ; rax <- rax * rax
    div     rcx
    mov     rax, rdx                ; rax <- rax mod p
    test    r10, r9
    jz      current_bit_not_set
    shl     rax, 1                  ; rax <- 2 * rax
    jc      mul2_overflow
    xor     rdx, rdx
    jmp     no_mul2_overflow
mul2_overflow:
    mov     rdx, 1
no_mul2_overflow:
    div     rcx
    mov     rax, rdx                ; rax <- rax mod p
current_bit_not_set:
    shr     r10, 1
    jmp     ble_loop
end_ble_loop:
    cmp     rax, 1
    je      ret_instruction         ; since rax == 1 we can leace immediately
    mov     r10, rcx
    dec     r10                     ; r10 <- p-1
sprp_loop:
    cmp     r10, rax
    je      return_true
    dec     r8
    jz      end_sprp_loop
    mul     rax                     ; rax <- rax * rax
    div     rcx
    mov     rax, rdx                ; rax <- rax mod p
    jmp     sprp_loop
end_sprp_loop:
    mov     al, 0
    jmp     ret_instruction
return_true:
    mov     al, 1
ret_instruction:
    ret
; ------------------------------------------------------------------------------------------------
is_2sprp endp
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
jacobi_symbol proc
; ------------------------------------------------------------------------------------------------
    test    rcx, rcx
    js      negative_a
    call    ujacobi_symbol
    ret
negative_a:
    neg     rcx
    mov     r10, rdx                ; backup m in r10 (it's OK because ujacobi_symbol doesn't use
                                    ; that register)
    call    ujacobi_symbol
    dec     r10
    shr     r10, 1
    test    r10, 1
    jz      do_not_neg_rax
    neg     rax
do_not_neg_rax:
    ret
    
; ------------------------------------------------------------------------------------------------
jacobi_symbol endp
; ------------------------------------------------------------------------------------------------


; ------------------------------------------------------------------------------------------------
; Computes the Jacobi symbol (a|m)
; for odd positive integer m and integer a >= 0
; ------------------------------------------------------------------------------------------------
; int64 ujacobi_symbol(uint64 a, uint64 m)
; a in rcx, m in rdx, return value in rax
; used registers: rax, rcx, rdx, r8, r9
; ------------------------------------------------------------------------------------------------
ujacobi_symbol proc
; ------------------------------------------------------------------------------------------------
    mov     r8, rdx                 ; backup m in r8
    xor     rdx, rdx
    mov     rax, rcx
    div     r8                      ; a mod m is in rdx
    mov     rcx, rdx                ; a <- a mod m
    mov     rdx, r8                 ; restore m
    mov     r8, 1                   ; t = 1
loop_a_not_zero:
    test    rcx, rcx
    jz      end_loop_a_not_zero
loop_a_even:
    test    rcx, 1
    jnz     end_loop_a_even
    shr     rcx, 1                  ; a <- a / 2
    mov     r9, rdx                 ; if m is 3 or 5 modulo 8 negate t
    and     r9, 7
    cmp     r9, 3
    je      neg_t
    cmp     r9, 5
    je      neg_t
    jmp     loop_a_even
neg_t:
    neg     r8
    jmp     loop_a_even
end_loop_a_even:
    xor     rcx, rdx                ;|
    xor     rdx, rcx                ;| -> swap a and m
    xor     rcx, rdx                ;|
    mov     r9, rcx                 ; if both a and m are 3 modulo 4 negate t
    and     r9, 3
    cmp     r9, 3
    jne     not_3_mod_4
    mov     r9, rdx
    and     r9, 3
    cmp     r9, 3
    jne     not_3_mod_4
    neg     r8
not_3_mod_4:
    mov     r9, rdx                 ; a <- a mod m
    xor     rdx, rdx
    mov     rax, rcx
    div     r9
    mov     rcx, rdx
    mov     rdx, r9
    jmp     loop_a_not_zero
end_loop_a_not_zero:
    cmp     rdx, 1
    je      return_t
    xor     rax, rax
    ret
return_t:
    mov     rax, r8
    ret 
; ------------------------------------------------------------------------------------------------
ujacobi_symbol endp
; ------------------------------------------------------------------------------------------------

end

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
