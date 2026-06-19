# Algoritmo Simplex — Implementação em Java (CLRS)

Implementação do Método Simplex em Java, baseada no Capítulo 29 do livro *Introduction to Algorithms* (CLRS) — Cormen, Leiserson, Rivest e Stein.

O algoritmo resolve problemas de Programação Linear na forma:

```
Maximizar   c^T · x
Sujeito a   A · x ≤ b   (m restrições)
                x ≥ 0   (n variáveis)
```

A implementação segue fielmente a notação do livro (Forma de Folgas, Regra de Bland para evitar ciclagem, e a operação PIVOT descrita na Seção 29.3).

## Autores

- Paulo Henrique Fonseca de Assis
- Rafael Franco

## Vídeo explicativo

Gravamos um vídeo explicando a teoria, a complexidade do algoritmo e um exemplo prático de execução passo a passo:

🎥 [Link do vídeo](https://drive.google.com/file/d/1Q_IVEGA74YFb6e0gey9I3oxA9AGvWVBP/view?usp=sharing)

## Referência

CORMEN, T. H.; LEISERSON, C. E.; RIVEST, R. L.; STEIN, C. *Introduction to Algorithms*, 3ª edição — Capítulo 29 (Linear Programming).
