#!/bin/bash

echo "Iniciando script de concatenação de arquivos FASTQ."
echo "Pressione Enter para confirmar que colou este script intencionalmente no terminal."
read

# Pergunta inicial
echo "Quantos pacientes você vai trabalhar? (ou digite 'exit' para sair)"
read PACIENTES
if [[ "$PACIENTES" == "exit" ]]; then exec bash; fi
if ! [[ "$PACIENTES" =~ ^[0-9]+$ ]]; then
  echo "Valor inválido. Digite um número inteiro."
  exec bash
fi

# Arquivos pareados?
echo "Os arquivos desses pacientes são pareados? (yes/no ou 'exit')"
read PAREADOS
if [[ "$PAREADOS" == "exit" ]]; then exec bash; fi

if [[ "$PAREADOS" != "yes" ]]; then
  echo "Este script só funciona com arquivos pareados. Encerrando."
  exec bash
fi

# Cálculo automático
TOTAL_ARQUIVOS=$((PACIENTES * 4))
R1_ARQUIVOS=$((PACIENTES * 2))
R2_ARQUIVOS=$((PACIENTES * 2))

echo ""
echo "Serão $TOTAL_ARQUIVOS arquivos no total:"
echo "- $R1_ARQUIVOS arquivos para R1"
echo "- $R2_ARQUIVOS arquivos para R2"
echo ""

# Caminho dos arquivos
echo "Informe o caminho COMPLETO onde estão os arquivos (sem barra ao final, ou 'exit' para sair):"
read CAMINHO
if [[ "$CAMINHO" == "exit" ]]; then exec bash; fi
if [[ ! -d "$CAMINHO" ]]; then
  echo "Diretório '$CAMINHO' não encontrado."
  exec bash
fi

echo "Arquivos encontrados em $CAMINHO:"
ls "$CAMINHO"
echo ""
echo "Todos os arquivos devem estar nesse diretório."
echo "Caso contrário, saia com 'exit' e mova os arquivos para um único local."
echo "Deseja continuar? (yes/no ou 'exit')"
read CONTINUE
if [[ "$CONTINUE" == "exit" ]]; then exec bash; fi
if [[ "$CONTINUE" != "yes" ]]; then echo "Cancelado."; exec bash; fi

# Arquivos R1
echo "Digite os nomes EXATOS dos $R1_ARQUIVOS arquivos de R1, separados por espaço:"
read -a R1_FILES
if [[ "${#R1_FILES[@]}" -ne "$R1_ARQUIVOS" ]]; then
  echo "Número incorreto de arquivos R1 fornecido."
  exec bash
fi

# Arquivos R2
echo "Digite os nomes EXATOS dos $R2_ARQUIVOS arquivos de R2, separados por espaço:"
read -a R2_FILES
if [[ "${#R2_FILES[@]}" -ne "$R2_ARQUIVOS" ]]; then
  echo "Número incorreto de arquivos R2 fornecido."
  exec bash
fi

# Pasta de saída
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$CAMINHO/1-arquivos_concatenados_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

echo ""
echo "Os arquivos concatenados serão salvos em: $OUTPUT_DIR"
echo "Deseja prosseguir com a concatenação? (yes/no ou 'exit')"
read CONFIRM_FINAL
if [[ "$CONFIRM_FINAL" == "exit" ]]; then exec bash; fi
if [[ "$CONFIRM_FINAL" != "yes" ]]; then echo "Cancelado."; exec bash; fi

# Concatenação
echo ""
echo "Iniciando concatenação..."
for ((i=0; i<PACIENTES; i++)); do
  R1_1=${R1_FILES[$((i*2))]}
  R1_2=${R1_FILES[$((i*2+1))]}
  R2_1=${R2_FILES[$((i*2))]}
  R2_2=${R2_FILES[$((i*2+1))]}

  PREFIX=$(echo "${R1_1}" | cut -d'_' -f1-2)

  echo "Concatenando arquivos da amostra: $PREFIX"
  cat "$CAMINHO/$R1_1" "$CAMINHO/$R1_2" > "$OUTPUT_DIR/${PREFIX}_R1_combined.fastq.gz"
  cat "$CAMINHO/$R2_1" "$CAMINHO/$R2_2" > "$OUTPUT_DIR/${PREFIX}_R2_combined.fastq.gz"

  echo "Concatenação concluída: $PREFIX"
done

echo ""
echo "Todos os arquivos foram concatenados com sucesso."
echo "Arquivos finais estão em: $OUTPUT_DIR"
exec bash
