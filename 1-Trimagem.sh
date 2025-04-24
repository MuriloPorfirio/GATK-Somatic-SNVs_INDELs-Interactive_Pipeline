#!/bin/bash

echo "⚠️  Script iniciado."
echo "Se você colou este script no terminal, pressione Enter para continuar (isso evita execução automática acidental)."
read

echo "Quantos pacientes você deseja processar? (ou digite 'exit' para sair)"
read N
if [[ "$N" == "exit" ]]; then exec bash; fi
if ! [[ "$N" =~ ^[0-9]+$ ]]; then
  echo "❌ Valor inválido. Digite apenas um número inteiro."
  exec bash
fi

TOTAL_FILES=$((N * 2))
echo "Isso corresponde a $TOTAL_FILES arquivos (R1 e R2 por paciente). Está correto? (yes/no ou 'exit')"
read CONFIRM_TOTAL
if [[ "$CONFIRM_TOTAL" == "exit" ]]; then exec bash; fi
if [[ "$CONFIRM_TOTAL" != "yes" ]]; then
  echo "Encerrando o script. Verifique as informações antes de prosseguir."
  exec bash
fi

echo "Todos os arquivos estão no mesmo diretório? (yes/no ou 'exit')"
read SAME_DIR
if [[ "$SAME_DIR" == "exit" ]]; then exec bash; fi
if [[ "$SAME_DIR" != "yes" ]]; then
  echo "Por favor, consolide todos os arquivos em um único diretório antes de continuar."
  exec bash
fi

echo "Informe o caminho completo do diretório com os arquivos (sem barra ao final, ou digite 'exit' para sair):"
read INPUT_DIR
if [[ "$INPUT_DIR" == "exit" ]]; then exec bash; fi
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Erro: diretório '$INPUT_DIR' não encontrado."
  exec bash
fi

echo "Listando arquivos .fastq ou .fastq.gz em '$INPUT_DIR'..."
find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \)

echo "Deseja processar todos os arquivos listados acima? (yes/no ou 'exit')"
read PROCESS_ALL
if [[ "$PROCESS_ALL" == "exit" ]]; then exec bash; fi

if [[ "$PROCESS_ALL" != "yes" ]]; then
  echo "Insira os nomes EXATOS dos arquivos que deseja trimar, separados por espaço (ou 'exit' para sair):"
  read -a FILES
  if [[ "${FILES[0]}" == "exit" ]]; then exec bash; fi
else
  FILES=($(find "$INPUT_DIR" -maxdepth 1 -type f \( -name "*.fastq" -o -name "*.fastq.gz" \) -printf "%f\n"))
fi

# Criação da pasta de saída
TIMESTAMP=$(date +"%d-%m-%Y_%Hh%Mm")
OUTPUT_DIR="$INPUT_DIR/arquivos_trimados_$TIMESTAMP"
mkdir -p "$OUTPUT_DIR"

echo "Os arquivos serão salvos em: $OUTPUT_DIR"
echo "Deseja continuar? (yes/no ou 'exit')"
read CONTINUE
if [[ "$CONTINUE" == "exit" ]]; then exec bash; fi
if [[ "$CONTINUE" != "yes" ]]; then
  echo "Processo cancelado."
  exec bash
fi

# Detectar recursos do sistema
TOTAL_CPUS=$(nproc)
SAFE_CPUS=$((TOTAL_CPUS / 2))
TOTAL_RAM_GB=$(free -g | awk '/^Mem:/ {print $2}')
DEFAULT_RAM_GB=4
SUGGESTED_RAM_GB=$((TOTAL_RAM_GB / 2))

echo ""
echo "O servidor possui $TOTAL_CPUS threads disponíveis."
echo "⚠️ Sugerimos até $SAFE_CPUS threads para evitar impacto em outros usuários."
echo "Quantas threads deseja usar? (ou 'exit')"
read THREADS
if [[ "$THREADS" == "exit" ]]; then exec bash; fi
if ! [[ "$THREADS" =~ ^[0-9]+$ ]]; then
  echo "Valor inválido para threads."
  exec bash
fi

echo ""
echo "Este script normalmente roda com ~${DEFAULT_RAM_GB}GB de RAM."
echo "🧠 O servidor possui $TOTAL_RAM_GB GB. Sugerimos até ${SUGGESTED_RAM_GB}GB."
echo "Quanto de RAM deseja alocar (GB)? (ou 'exit')"
read RAM
if [[ "$RAM" == "exit" ]]; then exec bash; fi
if ! [[ "$RAM" =~ ^[0-9]+$ ]]; then
  echo "Valor inválido para RAM."
  exec bash
fi

# Pergunta sobre uso de --user
echo ""
echo "Deseja rodar o Docker com um usuário específico (--user UID)?"
echo "Digite o UID desejado (ex: 1006), ou 'no' para continuar normalmente:"
read CUSTOM_UID
if [[ "$CUSTOM_UID" == "exit" ]]; then exec bash; fi
if [[ "$CUSTOM_UID" =~ ^[0-9]+$ ]]; then
  USE_USER="--user $CUSTOM_UID"
  echo "Docker será executado com: $USE_USER"
elif [[ "$CUSTOM_UID" == "no" ]]; then
  USE_USER=""
  echo "Docker será executado com usuário padrão (root dentro do container)."
else
  echo "❌ Entrada inválida. Digite um número de UID ou 'no'."
  exec bash
fi

# Resumo final
echo ""
echo "========= RESUMO FINAL ========="
echo "📁 Diretório de entrada: $INPUT_DIR"
echo "📂 Diretório de saída: $OUTPUT_DIR"
echo "📄 Arquivos a processar:"
for file in "${FILES[@]}"; do
  echo "- $file"
done
echo ""
echo "⚙️  Threads: $THREADS"
echo "🧠 RAM (informativo): $RAM GB"
echo "🔐 Docker UID: ${CUSTOM_UID}"
echo "================================"
echo "Podemos começar? (yes/no ou 'exit')"
read FINAL_CONFIRM
if [[ "$FINAL_CONFIRM" == "exit" ]]; then exec bash; fi
if [[ "$FINAL_CONFIRM" != "yes" ]]; then
  echo "Execução cancelada."
  exec bash
fi

# Início do processo
echo "Iniciando a trimagem com Docker..."

i=0
while [[ $i -lt ${#FILES[@]} ]]; do
  R1="${FILES[$i]}"
  R2="${FILES[$((i+1))]}"
  echo "Processando par:"
  echo "R1: $R1"
  echo "R2: $R2"

  docker run --rm \
    -v "$INPUT_DIR":/data \
    $USE_USER \
    biowardrobe2/trimgalore:v0.4.4 \
    trim_galore --paired "/data/$R1" "/data/$R2" -o "/data/arquivos_trimados_$TIMESTAMP"

  if [[ $? -ne 0 ]]; then
    echo "⚠️ Erro ao processar: $R1 e $R2"
  else
    echo "✅ Concluído: $R1 e $R2"
  fi

  i=$((i + 2))
done

echo "------------------------------------------------------"
echo "✅ Processo de trimagem finalizado."
echo "📁 Resultados em: $OUTPUT_DIR"
exec bash
