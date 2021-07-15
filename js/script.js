const form = document.getElementById('form');

form.addEventListener('submit', (e) => {
    e.preventDefault();
    let nome = document.getElementById('nome').value;
    let email = document.getElementById('email').value;
    let dados = {
        nome,
        email
    }
    let dadosConvertidos = JSON.stringify(dados);

    localStorage.setItem('lead', dadosConvertidos);

    let content = document.getElementById('content');
    let carregando = `<img src="assets/loading.gif" alt="loading" height="100" width="auto">`;
    let statusEnviado = `<h2>Cadastro realizado com sucesso!</h2>`;
    content.innerHTML = carregando;

    setTimeout(() => {
        content.innerHTML = statusEnviado;
    },1000);
    
})