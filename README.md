#  PDX Pipeline – Pré-processamento de FASTQ

Este repositório contém scripts em Bash para o pré-processamento de dados de sequenciamento genômico (FASTQ).

---

## ⚠️ Orientações importantes antes do uso

### Scripts interativos e longos

Todos os scripts deste repositório são **interativos**, com diversas perguntas ao usuário. Por isso, são visualmente grandes e **não devem ser colados diretamente em um terminal com screen**, pois:

- O conteúdo pode ser executado **instantaneamente**, sem tempo para resposta
- Isso pode travar ou fechar a sessão sem aviso

> **Recomendado:**
> 1. Crie um arquivo `.sh` dentro do screen usando o `nano`:
>    ```bash
>    nano meu_script.sh
>    ```
> 2. Cole o conteúdo do script no nano
> 3. Salve (`Ctrl+O`, `Enter`, `Ctrl+X`)
> 4. Torne o script executável:
>    ```bash
>    chmod +x meu_script.sh
>    ```
> 5. Execute dentro da screen com:
>    ```bash
>    ./meu_script.sh
>    ```

---

### Ordem correta: R1 sempre antes do R2

**Em todas as etapas do pipeline, até o alinhamento final**, é fundamental garantir que o arquivo **R1 venha antes de R2**, especialmente nos comandos:

- Trimagem (`trim_galore`)
- FastQC/MultiQC
- Alinhamento (Bowtie2, Parabricks, etc)

> Uma inversão nessa ordem pode afetar a qualidade da análise ou gerar erros silenciosos.

## Boas práticas

- Nomeie os arquivos com padrão claro (`_R1_`, `_R2_`, `_val_1`, `_val_2`)
- Execute os scripts **na ordem correta**
- Utilize `screen` ou `tmux` com cautela
- Revise os arquivos `.html` gerados
- Garanta sempre a ordem correta: **R1 primeiro, R2 depois**

---

## Autor
Murilo Porfírio de Aguiar e GPT
murilo.porfirio@yahoo.com
