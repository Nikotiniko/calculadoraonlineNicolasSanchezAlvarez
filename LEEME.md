# Web académica de Nicolás Sánchez Álvarez

Contenido: portada, 137 calculadoras (`calculadoras/`) y 8 herramientas interactivas (`herramientas/`). Es una web estática: no necesita servidor ni base de datos.

## Publicar en GitHub Pages

1. Descomprime el zip en tu ordenador. GitHub no descomprime archivos zip.
2. Abre tu repositorio `Nikotiniko/Myweb` en github.com.
3. Pulsa **Add file → Upload files**, selecciona **todo el contenido** de la carpeta (no la carpeta en sí) y arrástralo. Pulsa **Commit changes**. Si ya habías subido una versión anterior, los archivos se sustituyen.
4. Ve a **Settings → Pages**. En *Source* elige **Deploy from a branch**, rama **main** y carpeta **/(root)**, y pulsa **Save**. Con el plan gratuito, el repositorio tiene que ser **público**.
5. En unos minutos la web estará en:
   `https://nikotiniko.github.io/Myweb/`

Comprueba que se ha subido `.nojekyll`. Algunos sistemas ocultan los archivos que empiezan por punto. Si no aparece, créalo vacío con **Add file → Create new file**.

## Que Google la encuentre

- Da de alta la web en **Google Search Console** (gratuito) y envía el mapa del sitio: `sitemap.xml`. Es la forma más rápida de que Google indexe las 148 páginas.
- `robots.txt` solo lo leen los buscadores en la raíz de un dominio. En una dirección como la de arriba no tiene efecto, pero no estorba y servirá si más adelante usas un dominio propio o un repositorio `nikotiniko.github.io`.

## Si cambia la dirección de la web

La dirección completa aparece en `sitemap.xml`, en `robots.txt`, en los metadatos (`rel="canonical"` y `og:url`) de cada página y en el apartado «Cómo citar» de las herramientas. Si cambias el nombre del repositorio o usas un dominio propio, hay que actualizarla en esos sitios. Pídemelo y lo regenero todo de una vez.

## Pendiente de completar en la portada

- Publicaciones con DOI y enlaces a ORCID, Google Scholar y ResearchGate.
- Manuales de metodología y enlaces de Amazon de cada libro.
- Proyectos financiados (entidad, referencia y periodo).
