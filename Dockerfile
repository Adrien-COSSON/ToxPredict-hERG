# Image de base Python 3.10
FROM python:3.10

# ── SPECIFIC TO HUGGING FACE ─────────────────────────────────────────
# Create a user "user" with UID 1000 (required by HF Spaces)
RUN useradd -m -u 1000 user

# Pass from "root" to "user"
USER user

# Définir HOME et ajouter ~/.local/bin au PATH
# Nécessaire pour que pip install --user soit accessible
ENV HOME=/home/user \
    PATH=/home/user/.local/bin:$PATH
# ────────────────────────────────────────────────────────────────────

# Working directory = /home/user/app
WORKDIR $HOME/app

# Copier tous les fichiers locaux avec "user" comme propriétaire
# --chown=user est indispensable sur HF Spaces pour éviter les erreurs de permission
COPY --chown=user . $HOME/app

# Installing dependencies
RUN pip install -r requirements.txt

# Run MLflow server
# $PORT is provide by Huggin Face Spaces (7860 par défaut)
# Avantage de la variable : le code est plus portable si on change d'infra
# Modifier '$PORT' par le numero du port
CMD mlflow server -p $PORT \
    --host 0.0.0.0 \
    --allowed-hosts "adrien-cosson-toxpredict-mlflow.hf.space" \
    --cors-allowed-origins "https://adrien-cosson-toxpredict-mlflow.static.hf.space" \
    --backend-store-uri $BACKEND_STORE_URI \
    --default-artifact-root $ARTIFACT_STORE_URI