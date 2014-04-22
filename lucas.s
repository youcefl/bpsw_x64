; ------------------------------------------------------------------------------------------------
; Creation date: 2014.04.21
; Creator: Youcef Lemsafer
; Authors: Youcef Lemsafer
; ------------------------------------------------------------------------------------------------

; ------------------------------------------------------------------------------------------------
bits 64
default rel
; ------------------------------------------------------------------------------------------------

; ------------------------------------------------------------------------------------------------
SECTION .code
; ------------------------------------------------------------------------------------------------

global is_slprp

; ------------------------------------------------------------------------------------------------
; We assume that n is odd > 1 and that Q = (1 - D)/4
; where D is the first element in {5, -7, 9, -11, ...}
; such that the Jacobi symbol (D/n) = -1.
; ------------------------------------------------------------------------------------------------
; bool is_slprp(uint64 n, int64 D);
; ------------------------------------------------------------------------------------------------
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
is_slprp:
    sub     rsp, 28h
    mov     rcx, rdi    ; this first two lines is to make
    mov     rdx, rsi    ; it look like Windows, they have to be removed
                        ; if we are on Windows
    mov     [rsp], rcx
    mov     [rsp + 8h], rdx
    mov     rax, rdx
    neg     rax
    add     rax, 1
    sar     rax, 2      ; after that Q is in rax
    mov     [rsp + 10h], rax

    ; Factor n + 1 into d.2^s, d in r8 and s in r9
    mov     r10, 8000000000000000h
    mov     r8, rcx     ; initialize d = n + 1
    add     r8, 1       ; we do not test for overflow here since 2^64-1 =
                        ; 3.5.17.257.641.65537.P7 cannot end here (because
                        ; we assume that trial division by small primes has
                        ; been performed).
    mov     r9, 0       ; initialize s = 0
div2:
    test    r8, 1
    jnz     div2_end
    shr     r8, 1
    shr     r10, 1
    inc     r9
    jmp     div2
div2_end:               ; at this point d is in r8 and s is in r9

    ; Compute U(d) mod n and V(d) mod n

adj_d_msbmask:
    test    r8, r10
    jnz     adj_d_msbmask_end
    shr     r10, 1
    jmp     adj_d_msbmask
adj_d_msbmask_end:

    ; initialization
    mov     rcx, 0      ; U = U(0) = 0
    mov     r11, 2      ; V = V(0) = 2
    mov     qword [rsp + 18h], 1    ; q is the k-th power of Q
loop_d_msbmask:         ; loop until r10 is 0
    test    r10, r10
    jz      loop_d_msbmask_end
    mov     rax, rcx    ; compute UV
    imul     r11
    idiv    qword [rsp]
    mov     rcx, rdx    ; U = UV mod n

    mov     rax, r11    ; compute V^2 - 2Q^k = V^2 - 2q (mod n)
    imul     r11
    idiv     qword [rsp]
    mov     r11, rdx    ; r11 <- V^2 (mod n)
    mov     rax, 2
    imul     qword [rsp + 18h]
    idiv    qword [rsp]
    sub     r11, rdx    ; r11 <- r11 - 2q (mod n)

    mov     rax, qword [rsp + 18h]      ; update q for next round:
    imul     rax
    idiv    qword [rsp]
    mov     qword [rsp + 18h], rdx      ; q <- q^2 (mod n)

    test    r10, r8
    jz      current_bit_not_set

    ; U(2k+1) = (U(2k) + V(2k)) / 2
    mov     rax, rcx
    add     rax, r11                    ; @todo: what about possible overflow ?
    test    rax, 1
    jz      uv_sum_is_even
    add     rax, qword [rsp]
uv_sum_is_even:
    sar     rax, 1
    mov     qword [rsp + 20h], rax      ; backup U(2k+1), we still need U(2k)

    ; V(2k+1) = (D.U(2k) + V(2k)) / 2
    mov     rax, rcx                    ; rax <- U(2k)
    imul    qword [rsp + 8h]            ; rax <- U(2k) * D
    add     rax, r11                    ; @todo: what about possible overflow ?
    test    rax, 1
    jz      duv_sum_is_even
    add     rax, qword [rsp]
duv_sum_is_even:
    sar     rax, 1                      ; rax <- rax / 2
    mov     rcx, qword [rsp + 20h]      ; rcx <- U(2k+1)
    mov     r11, rax                    ; r11 <- V(2k+1)

    ; q <- q * Q (mod n)
    mov     rax, qword [rsp + 18h]      ; rax <- q
    imul    qword [rsp + 10h]           ; multiply by Q
    idiv    qword [rsp]
    mov     qword [rsp + 18h], rdx      ; q <- qword [rsp + 18h]

current_bit_not_set:
    shr     r10, 1
    jmp     loop_d_msbmask

loop_d_msbmask_end:

    ; Now U(d) is in rcx and V(d) is in r11
    test    rcx, rcx
    jz      is_slpsp_true

    xor     rax, rax
    jmp     is_slpsp_end

is_slpsp_true:
    mov     rax, 1

is_slpsp_end:
    add     rsp, 28h
    ret




