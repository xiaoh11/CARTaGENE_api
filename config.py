import os.path

# Required configuration
#MONGO_URI = 'mongodb://localhost:27017/test'
MONGO_URI = 'mongodb://127.0.0.1:27017/test'
BRAVO_API_PAGE_LIMIT = 10000

BASE_DIR = os.path.join(os.sep, 'home', 'xiaoh11', "projects", "rrg-vmooser", "xiaoh11", 'bravo', 'data', 'runtime')
# BASE_DIR = os.path.join(os.sep, 'home', 'xiaoh', 'bravo', 'data', 'runtime')
COVERAGE_DIR = os.path.join(BASE_DIR, 'coverage')
SEQUENCES_DIR = os.path.join(BASE_DIR, 'crams')
SEQUENCES_CACHE_DIR = os.path.join(BASE_DIR, 'cache')
REFERENCE_SEQUENCE = os.path.join(BASE_DIR, 'reference', 'hs38DH.fa')
CLINVAR_DIR = os.path.join(BASE_DIR, 'clinvar')
CLINVAR_VCF = os.path.join(BASE_DIR, 'clinvar', 'clinvar_20231007.vcf.gz')

# Optional configuration
LOGIN_DISABLED = True
SESSION_SECRET = b'deadbeef0123456789'
CORS_ORIGINS = ['http://localhost:8083']

# Config for using Google OAuth
GOOGLE_CLIENT_ID = "your google oauth client id"
GOOGLE_CLIENT_SECRET = "your google oauth client secret"
GOOGLE_DISCOVERY_URL = "https://accounts.google.com/.well-known/openid-configuration"
