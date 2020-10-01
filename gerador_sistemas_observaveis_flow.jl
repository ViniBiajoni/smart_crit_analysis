clearconsole()
# Caminho = "D:/Programas Julia/Gerador Sistemas Observaveis"
 Caminho = "C:/Users/vinib/OneDrive/Documentos/Doutorado/Tec. Inteligentes em Sist Pot/Monta_Sist_Observaveis"
# Caminho = "C:/Users/Rafael/OneDrive/Julia/meas-scheme-gen"

cd(Caminho)

using LinearAlgebra
using Dates
using DelimitedFiles


# Monta vetor de medidas completo
function vetMed(A)

    m = sum(UpperTriangular(A)) # Número de ramos
    n = size(A, 1)
    nMedMax = n + 2 * m     # Número máximo de medidas possível
    medidas = zeros(Int32, nMedMax, 3) # Vetor que organiza medidores.

    global med = 0

    for de = 1:n
        for para = 1:n
            if (A[de, para] == 1)
                med = med + 1
                medidas[med, 1] = Int(de)
                medidas[med, 2] = Int(para)
            end
        end
    end #flUXO

    numflow = med

    for i = 1:n
        med = med + 1
        medidas[med, 1] = Int(i)
        medidas[med, 2] = Int(i)
    end #INJEÇÃO


    return (medidas, numflow)
end


# Calcula matriz jacobiano
function jacobiana(medidas, A, ref)
    nMed = Int(sum(medidas[:, 3]))
    n = size(A, 1)
    H = zeros(nMed + 1, n)
    H[nMed+1, ref] = 1
    ind = 1
    med = 1
    while med <= nMed
        if medidas[ind, 3] == 1
            de = Int(medidas[ind, 1])
            para = Int(medidas[ind, 2])
            if de != para
                for l = 1:n
                    if (l == de)
                        H[med, l] = 1
                    elseif (l == para)
                        H[med, l] = -1
                    end
                end
            else
                nbc = Int(sum(A[de, :]))
                for l = 1:n
                    if (l == de)
                        H[med, l] = nbc
                    else
                        H[med, l] = -1 * (A[de, l])
                    end
                end
            end
            med = med + 1
        end
        ind = ind + 1
    end

    return H
end

## Selecion a Rede Desejada
#net_connection = input("Tamanho da Rede")
net_connection = 14

if (net_connection == 14)

    A = readdlm("Conexao14Bus.txt", ',', Int, '\n')

end

if (net_connection == 30)

    A = readdlm("Conexao30Bus.txt", ',', Int, '\n')

end

if (net_connection == 118)

    A = readdlm("Conexao118Bus.txt", ',', Int, '\n')

end


function  Gera_Cov(A,nmed)
    (full_med, numflow) = vetMed(A) # monta a versão
    # Constroi o vetor que organiza as medidas disponíveis.
    # Da forma: [|barra de| |barra para| |ligado ou desligado|]
    # Ativa, na matriz medidas, as medidas em medidasAtivas
    #(medidas,actIdx)=altMed(medidas,medidasAtivas,1)
    #println(numflow)

    n_max_med = size(full_med,1)
    H=jacobiana(full_med,A,1)
    G = transpose(H) * H
    med=0
    sorteio = collect(1:1:numflow)

    while  (abs(det(G)) < 1E-4 || med<nmed) && numflow != 0


        #Sorteia
        val_sorteado = rand(1:numflow)
        #Armazena a medida

        full_med[sorteio[val_sorteado],3] = 1 #liguei a medida de fluxo
        #Deleta Medida ja Sorteada
        sorteio[val_sorteado]=0
        #Atualiza
        sorteio = filter!(!iszero, sorteio)

        #display(sorteio)
        numflow = numflow-1
        H=jacobiana(full_med,A,1)
        G = transpose(H) * H
        med=med+1
        #Q,R = qr(G)

    end
    println(med)
        #Gera a matriz de Covariancias
        E= I - H*inv(G)*transpose(H)

        return(E,full_med)
end


numero_casos = 1 # numero de casos a serem gerados

#Geracao de multiplos casos
for i=1:numero_casos

    ### Qtde minima de medidas para cada plano a ser Gerador
    nmed= 26 #14 Barras
    #nmed= 55 #30 Barras
    #nmed= 225 #118 Barras
    (E,full_med)=Gera_Cov(A,nemd)
    nmedidores=sum(full_med[:,3])

    ###Imprime a matriz de Covariâncias
    date=string(Date(now()),"_",hour(now()),"h",minute(now()),"m")
    nomeE=string("E",size(A,1),"b",nmedidores,"m",date,".txt")

    open(nomeE, "w") do io
        writedlm(io, E, ',')
    end

    ####Imprime o arquivo de medidas

    nome_sistMed = string("Med_Plan",size(A,1),"b",nmedidores,"m",date,".txt")
    open(nome_sistMed, "w") do io
        writedlm(io, full_med, ',')
    end

end
